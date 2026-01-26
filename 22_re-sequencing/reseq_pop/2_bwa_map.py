#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 FASTQ 测序数据比对到参考基因组，并进行标准的 BAM 处理。

本脚本执行以下步骤：
1. BWA MEM: 将双端测序数据比对到参考基因组
2. Samtools sort: 对 BAM 文件按坐标排序
3. Samtools/GATK MarkDuplicates: 标记 PCR 重复（提高变异检测准确性）
4. Samtools index: 为最终 BAM 文件建立索引

软件依赖:
- bwa
- samtools
- gatk (可选，用于更精确的 MarkDuplicates)
"""

import os
import argparse
import subprocess
import sys
import glob
import concurrent.futures


def run_command(command, description, sample_name, exit_on_fail=True):
    """运行命令并打印状态，带有样本名标识"""
    prefix = f"[{sample_name}]"
    print(f"\n{prefix} [子步骤] {description}")
    print(f"{prefix} [命令] {command}")
    
    # 捕获输出以防日志混乱，或者直接打印（并行时会交错但带有前缀）
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print(f"{prefix} [错误] 命令执行失败，返回码: {result.returncode}")
        if exit_on_fail:
            # 在多进程中直接退出可能不会停止整个脚本，但会报异常
            raise RuntimeError(f"Sample {sample_name} failed at step: {description}")
        return False
    return True


def find_samples(folder, plus_tag, minus_tag):
    """在指定文件夹中查找所有样本"""
    samples = []
    for filepath in glob.glob(os.path.join(folder, f"*{plus_tag}")):
        filename = os.path.basename(filepath)
        sample_name = filename.replace(plus_tag, "")
        
        # 验证对应的反向读段文件是否存在
        minus_file = os.path.join(folder, f"{sample_name}{minus_tag}")
        if os.path.exists(minus_file):
            samples.append({
                'name': sample_name,
                'r1': filepath,
                'r2': minus_file
            })
        else:
            print(f"警告: 样本 '{sample_name}' 缺少反向读段文件 ({minus_tag})，跳过。")
    
    return samples


def process_sample(sample, reference, threads, output_folder, use_gatk_markdup):
    """处理单个样本的完整比对流程"""
    
    sample_name = sample['name']
    r1_path = sample['r1']
    r2_path = sample['r2']
    
    # 定义进度输出前缀
    p = f"[{sample_name}]"
    
    print(f"\n{p} {'='*40}")
    print(f"{p} 开始处理: R1={os.path.basename(r1_path)}, R2={os.path.basename(r2_path)}")
    print(f"{p} {'='*40}")
    
    # 定义中间文件和输出文件路径
    raw_bam = os.path.join(output_folder, f"{sample_name}.raw.bam")
    sorted_bam = os.path.join(output_folder, f"{sample_name}.sorted.bam")
    final_bam = os.path.join(output_folder, f"{sample_name}.bam")
    markdup_metrics = os.path.join(output_folder, f"{sample_name}.markdup_metrics.txt")
    
    try:
        # 步骤 1: BWA MEM 比对
        rg_string = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tLB:WGS\\tPL:ILLUMINA"
        bwa_cmd = (
            f"bwa mem -t {threads} -M "
            f"-R '{rg_string}' "
            f"{reference} {r1_path} {r2_path} | "
            f"samtools view -bS - > {raw_bam}"
        )
        run_command(bwa_cmd, "BWA MEM 比对", sample_name)
        
        # 步骤 2: Samtools sort 排序
        sort_cmd = f"samtools sort -@ {threads} -o {sorted_bam} {raw_bam}"
        run_command(sort_cmd, "Samtools 排序", sample_name)
        
        # 删除原始 BAM
        if os.path.exists(raw_bam):
            os.remove(raw_bam)
        
        # 步骤 3: 标记 PCR 重复
        if use_gatk_markdup:
            markdup_cmd = (
                f"gatk MarkDuplicates "
                f"-I {sorted_bam} "
                f"-O {final_bam} "
                f"-M {markdup_metrics} "
                f"--REMOVE_DUPLICATES false "
                f"--CREATE_INDEX true"
            )
            success = run_command(markdup_cmd, "GATK MarkDuplicates (标记PCR重复)", sample_name, exit_on_fail=False)
            
            if not success:
                print(f"{p} [信息] GATK MarkDuplicates 失败，回推到 samtools markdup...")
                use_gatk_markdup = False
        
        if not use_gatk_markdup:
            markdup_cmd = f"samtools markdup -@ {threads} {sorted_bam} {final_bam}"
            run_command(markdup_cmd, "Samtools markdup (标记PCR重复)", sample_name)
            
            index_cmd = f"samtools index {final_bam}"
            run_command(index_cmd, "Samtools 索引", sample_name)
        
        # 删除排序后的中间 BAM
        if os.path.exists(sorted_bam):
            os.remove(sorted_bam)
        
        print(f"\n{p} [完成] 处理完毕。输出: {final_bam}")
        return final_bam
    except Exception as e:
        print(f"\n{p} [严重错误] 样本处理中断: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="将 FASTQ 数据比对到参考基因组并进行并行 BAM 处理。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    required = parser.add_argument_group("必需参数")
    required.add_argument('-r', '--reference', required=True, help="参考基因组 FASTA 文件")
    required.add_argument('-p', '--plus-tag', required=True, help="正向读段后缀 (如 _1.fq.gz)")
    required.add_argument('-m', '--minus-tag', required=True, help="反向读段后缀 (如 _2.fq.gz)")
    
    optional = parser.add_argument_group("可选参数")
    optional.add_argument('-t', '--threads', type=int, default=4, help="每个样本使用的线程数 (默认: 4)")
    optional.add_argument('-j', '--parallel-samples', type=int, default=1, help="同时处理的样本数 (默认: 1)")
    optional.add_argument('-f', '--folder', default=".", help="FASTQ 文件夹")
    optional.add_argument('-o', '--output', default=None, help="BAM 输出文件夹")
    optional.add_argument('--use-gatk-markdup', action='store_true', help="使用 GATK MarkDuplicates")
    
    args = parser.parse_args()
    
    input_folder = os.path.abspath(args.folder)
    output_folder = os.path.abspath(args.output) if args.output else input_folder
    reference = os.path.abspath(args.reference)
    
    # 检查环境和文件... (省略部分检查)
    os.makedirs(output_folder, exist_ok=True)
    samples = find_samples(input_folder, args.plus_tag, args.minus_tag)
    
    if not samples:
        print("未找到有效样本，退出。")
        sys.exit(1)

    # 核心并行逻辑
    num_parallel = min(args.parallel_samples, len(samples))
    threads_per_sample = args.threads
    total_threads = num_parallel * threads_per_sample
    
    print(f"\n[任务启动] 样本总数: {len(samples)}")
    print(f"[资源配置] 并行样本数: {num_parallel}")
    print(f"[资源配置] 每样本线程: {threads_per_sample}")
    print(f"[资源配置] 预计峰值总线程: {total_threads}")
    
    output_bams = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel) as executor:
        # 提交所有任务
        future_to_sample = {
            executor.submit(
                process_sample, 
                sample, 
                reference, 
                threads_per_sample, 
                output_folder, 
                args.use_gatk_markdup
            ): sample['name'] for sample in samples
        }
        
        # 收集结果
        for future in concurrent.futures.as_completed(future_to_sample):
            sample_name = future_to_sample[future]
            try:
                bam_path = future.result()
                if bam_path:
                    output_bams.append(bam_path)
            except Exception as exc:
                print(f"{sample_name} 生成了未预料的异常: {exc}")

    print("\n" + "="*60)
    print(f"所有任务处理尝试完成！成功: {len(output_bams)}/{len(samples)}")
    print("="*60)
    for bam in output_bams:
        print(f"  - {bam}")


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()