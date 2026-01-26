#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用 GATK HaplotypeCaller 进行变异检测（并行优化版）。

优化策略：
- Scatter-Gather: 将基因组按染色体拆分 (Scatter)，并行运行 HaplotypeCaller。
- 随后将每个样本的分片结果合并 (Gather)。
- 这种方式可以防止单条大染色体阻塞整个流程，并显著提高在大算力服务器上的利用率。

本脚本执行以下步骤：
1. 读取参考基因组索引 (.fai)，获取染色体列表
2. 对每个样本的每个染色体生成任务，并行运行 GATK HaplotypeCaller
3. 使用 GATK GatherVcfs 合并每个样本的分片 gVCF/VCF
4. (可选) 使用 CombineGVCFs + GenotypeGVCFs 合并所有样本

软件依赖:
- gatk (版本 4.x)
- python3
"""

import os
import argparse
import subprocess
import sys
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_command(command, description, exit_on_fail=True):
    """运行命令并打印状态"""
    # 简略日志以免刷屏
    # print(f"  [DEBUG] {command}") 
    
    result = subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"\n  [错误] {description} 失败")
        print(f"  [具体命令] {command}")
        print(f"  [Stderr] {result.stderr.decode('utf-8')}")
        if exit_on_fail:
            sys.exit(1)
        return False
    return True

def get_chromosomes(fai_path):
    """从 .fai 文件获取染色体列表"""
    chroms = []
    if not os.path.exists(fai_path):
        print(f"错误: 找不到 .fai 索引文件: {fai_path}")
        sys.exit(1)
        
    with open(fai_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                chroms.append(parts[0])
    return chroms

def call_variant_shard(bam_path, reference, output_folder, chrom, emit_gvcf, gatk_threads):
    """对单个样本的单条染色体进行变异检测"""
    sample_name = os.path.basename(bam_path).replace('.bam', '')
    
    # 定义分片输出文件名
    # 例如: vcf_output/shards/sample1/sample1.chr1.g.vcf.gz
    shard_dir = os.path.join(output_folder, "shards", sample_name)
    os.makedirs(shard_dir, exist_ok=True)
    
    suffix = ".g.vcf.gz" if emit_gvcf else ".vcf.gz"
    output_file = os.path.join(shard_dir, f"{sample_name}.{chrom}{suffix}")
    
    # 如果文件已存在且不为空，跳过（断点续传简单实现）
    if os.path.exists(output_file) and os.path.getsize(output_file) > 1000:
        return output_file

    erc_option = "--emit-ref-confidence GVCF" if emit_gvcf else ""
    
    # GATK 命令
    # -L: 指定染色体/区间
    cmd = (
        f"gatk HaplotypeCaller "
        f"-R {reference} "
        f"-I {bam_path} "
        f"-O {output_file} "
        f"-L {chrom} "
        f"{erc_option} "
        f"--native-pair-hmm-threads {gatk_threads}"
    )
    
    desc = f"{sample_name} - {chrom}"
    run_command(cmd, desc)
    
    # print(f"  [完成] {desc}")
    return output_file

def gather_shards(sample_name, shards, output_folder, emit_gvcf, reference):
    """合并单个样本的所有染色体分片"""
    print(f"\n[合并] 正在合并样本 {sample_name} 的 {len(shards)} 个分片...")
    
    suffix = ".g.vcf.gz" if emit_gvcf else ".vcf.gz"
    final_output = os.path.join(output_folder, f"{sample_name}{suffix}")
    
    # 构建输入参数 -I shard1 -I shard2 ...
    # 确保分片顺序正确（虽然 GatherVcfs 通常按基因组顺序，但传入有序列表更安全）
    # shards 列表在传入前应已排序
    inputs = " ".join([f"-I {s}" for s in shards])
    
    cmd = (
        f"gatk GatherVcfsCloud "
        f"{inputs} "
        f"-O {final_output} "
        f"--ignore-safety-checks" # 忽略可能因索引差异导致的检查
    )
    
    # GATK 4 中 GatherVcfsCloud 比较通用，GatherVcfs 也可以
    # 如果 GatherVcfsCloud 不存在，回退到 GatherVcfs
    if subprocess.run("gatk GatherVcfsCloud --help", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode != 0:
         cmd = cmd.replace("GatherVcfsCloud", "GatherVcfs")

    run_command(cmd, f"合并样本 {sample_name}")
    
    # 为合并后的文件建立索引
    run_command(f"gatk IndexFeatureFile -I {final_output}", f"建立索引 {sample_name}")
    
    return final_output

def combine_gvcfs(gvcf_files, reference, output_folder):
    """(联合调用模式) 合并所有样本的 gVCF"""
    print(f"\n{'='*60}")
    print("[步骤] 合并所有样本的 gVCF 文件 (CombineGVCFs)")
    print('='*60)
    
    combined_gvcf = os.path.join(output_folder, "cohort.g.vcf.gz")
    
    # 写入 list 文件以避免命令行过长
    list_file = os.path.join(output_folder, "gvcfs.list")
    with open(list_file, 'w') as f:
        for gvcf in gvcf_files:
            f.write(gvcf + "\n")
            
    combine_cmd = (
        f"gatk CombineGVCFs "
        f"-R {reference} "
        f"-V {list_file} "
        f"-O {combined_gvcf}"
    )
    
    run_command(combine_cmd, "GATK CombineGVCFs")
    return combined_gvcf

def genotype_gvcfs(combined_gvcf, reference, output_folder):
    """对合并的 gVCF 进行基因型分析"""
    print(f"\n{'='*60}")
    print("[步骤] 进行联合基因型分析 (Joint Genotyping)")
    print('='*60)
    
    final_vcf = os.path.join(output_folder, "all_samples.vcf.gz")
    
    genotype_cmd = (
        f"gatk GenotypeGVCFs "
        f"-R {reference} "
        f"-V {combined_gvcf} "
        f"-O {final_vcf} "
        f"--include-non-variant-sites"
    )
    
    run_command(genotype_cmd, "GATK GenotypeGVCFs")
    return final_vcf

def main():
    parser = argparse.ArgumentParser(
        description="GATK HaplotypeCaller 优化版 (Scatter-Gather)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    required = parser.add_argument_group("必需参数")
    required.add_argument('-r', '--reference', required=True, help="参考基因组 FASTA")
    
    optional = parser.add_argument_group("可选参数")
    optional.add_argument('-b', '--bam-folder', default=".", help="BAM 目录")
    optional.add_argument('-o', '--output', default="./vcf_output", help="输出目录")
    optional.add_argument('--joint-calling', action='store_true', help="联合变异检测流程")
    optional.add_argument('--single-vcf', action='store_true', help="仅生成单样本 VCF")
    
    optional.add_argument('-t', '--threads', type=int, default=20, 
                          help="总并发任务数 (建议设置较高，如 20-50)")
    optional.add_argument('--gatk-threads', type=int, default=1, 
                          help="每个 GATK 进程的线程数 (建议 1-2，把资源留给并发任务)")
    
    args = parser.parse_args()
    
    # 路径处理
    bam_folder = os.path.abspath(args.bam_folder)
    output_folder = os.path.abspath(args.output)
    reference = os.path.abspath(args.reference)
    
    # 检查参考文件
    if not os.path.exists(reference):
        print(f"错误: 参考文件不存在: {reference}")
        sys.exit(1)
    
    fai_file = reference + ".fai"
    if not os.path.exists(fai_file):
        print(f"错误: 缺少 .fai 索引，请先运行 samtools faidx")
        sys.exit(1)
        
    dict_file = reference.rsplit('.', 1)[0] + ".dict"
    if not os.path.exists(dict_file):
        print(f"错误: 缺少 .dict 文件")
        sys.exit(1)

    # 获取染色体
    chroms = get_chromosomes(fai_file)
    print(f"\n[信息] 检测到 {len(chroms)} 条染色体/Contigs，将进行并行分片处理。")
    
    os.makedirs(output_folder, exist_ok=True)
    
    # 获取 BAM
    bam_files = sorted(glob.glob(os.path.join(bam_folder, "*.bam")))
    bam_files = [b for b in bam_files if not any(x in b for x in ['.raw.', '.sorted.'])]
    
    if not bam_files:
        print(f"错误: 未找到 BAM 文件")
        sys.exit(1)
        
    print(f"[信息] 待处理样本数: {len(bam_files)}")
    
    emit_gvcf = args.joint_calling or (not args.single_vcf)
    
    # ==========================
    # Phase 1: Scatter (Mapping)
    # ==========================
    print(f"\n" + "="*60)
    print(f"[阶段 1] 并行变异检测 (Scatter)")
    print(f"  - 总并发限制: {args.threads}")
    print(f"  - 单任务线程: {args.gatk_threads}")
    print(f"  - 总任务数: {len(bam_files) * len(chroms)}")
    print("="*60)
    
    # 准备所有任务
    all_tasks = []
    # 结构: (sample_index, chrom_index, bam_path, chrom)
    # 用于排序结果
    
    for bam in bam_files:
        for chrom in chroms:
            all_tasks.append((bam, chrom))
            
    shard_results = {} # {bam_path: {chrom: output_file}}
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        # submit all
        future_to_task = {
            executor.submit(
                call_variant_shard, 
                bam, reference, output_folder, chrom, emit_gvcf, args.gatk_threads
            ): (bam, chrom)
            for bam, chrom in all_tasks
        }
        
        completed_count = 0
        total_count = len(all_tasks)
        
        for future in as_completed(future_to_task):
            bam, chrom = future_to_task[future]
            try:
                out_file = future.result()
                if bam not in shard_results:
                    shard_results[bam] = {}
                shard_results[bam][chrom] = out_file
                
                completed_count += 1
                if completed_count % 10 == 0:
                    print(f"  [进度] {completed_count}/{total_count} 分片完成...")
                    
            except Exception as e:
                print(f"[错误] 任务失败 ({os.path.basename(bam)} - {chrom}): {e}")
                sys.exit(1)

    print(f"\n[信息] 所有分片检测完成。")

    # ==========================
    # Phase 2: Gather (Merging)
    # ==========================
    print(f"\n" + "="*60)
    print(f"[阶段 2] 合并分片 (Gather)")
    print("="*60)
    
    final_sample_files = []
    
    for bam in bam_files:
        sample_name = os.path.basename(bam).replace('.bam', '')
        # 按染色体顺序收集分片文件
        sample_shards = []
        missing = False
        for chrom in chroms:
            if chrom in shard_results.get(bam, {}):
                sample_shards.append(shard_results[bam][chrom])
            else:
                print(f"[错误] 样本 {sample_name} 缺少染色体 {chrom} 的结果！")
                missing = True
        
        if missing:
            sys.exit(1)
            
        merged_file = gather_shards(sample_name, sample_shards, output_folder, emit_gvcf, reference)
        final_sample_files.append(merged_file)

    # ==========================
    # Phase 3: Joint Calling
    # ==========================
    if args.joint_calling and len(final_sample_files) > 1:
        combine_gvcfs(final_sample_files, reference, output_folder)
        final_vcf = genotype_gvcfs(os.path.join(output_folder, "cohort.g.vcf.gz"), reference, output_folder)
        print(f"\n[完成] 联合 VCF: {final_vcf}")
    else:
        print(f"\n[完成] 单样本处理完毕。")

if __name__ == "__main__":
    main()
