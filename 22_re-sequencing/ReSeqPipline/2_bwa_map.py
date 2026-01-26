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


def run_command(command, description, exit_on_fail=True):
    """运行命令并打印状态"""
    print(f"\n  [子步骤] {description}")
    print(f"  [命令] {command}")
    
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print(f"  [错误] 命令执行失败，返回码: {result.returncode}")
        if exit_on_fail:
            sys.exit(1)
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
    
    print(f"\n{'='*60}")
    print(f"[样本] {sample_name}")
    print(f"  R1: {os.path.basename(r1_path)}")
    print(f"  R2: {os.path.basename(r2_path)}")
    print('='*60)
    
    # 定义中间文件和输出文件路径
    raw_bam = os.path.join(output_folder, f"{sample_name}.raw.bam")
    sorted_bam = os.path.join(output_folder, f"{sample_name}.sorted.bam")
    final_bam = os.path.join(output_folder, f"{sample_name}.bam")
    markdup_metrics = os.path.join(output_folder, f"{sample_name}.markdup_metrics.txt")
    
    # 步骤 1: BWA MEM 比对
    # -M: 将短的 split hits 标记为 secondary（兼容 Picard）
    # -R: 添加 Read Group 信息（GATK 必须）
    rg_string = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tLB:WGS\\tPL:ILLUMINA"
    bwa_cmd = (
        f"bwa mem -t {threads} -M "
        f"-R '{rg_string}' "
        f"{reference} {r1_path} {r2_path} | "
        f"samtools view -bS - > {raw_bam}"
    )
    run_command(bwa_cmd, "BWA MEM 比对")
    
    # 步骤 2: Samtools sort 排序
    sort_cmd = f"samtools sort -@ {threads} -o {sorted_bam} {raw_bam}"
    run_command(sort_cmd, "Samtools 排序")
    
    # 删除原始 BAM
    os.remove(raw_bam)
    
    # 步骤 3: 标记 PCR 重复
    if use_gatk_markdup:
        # 使用 GATK MarkDuplicates（更精确，内存消耗较大）
        markdup_cmd = (
            f"gatk MarkDuplicates "
            f"-I {sorted_bam} "
            f"-O {final_bam} "
            f"-M {markdup_metrics} "
            f"--REMOVE_DUPLICATES false "
            f"--CREATE_INDEX true"
        )
        success = run_command(markdup_cmd, "GATK MarkDuplicates (标记PCR重复)", exit_on_fail=False)
        
        if not success:
            print("  [信息] GATK MarkDuplicates 失败，回退到 samtools markdup...")
            use_gatk_markdup = False
    
    if not use_gatk_markdup:
        # 使用 samtools markdup（轻量级替代方案）
        markdup_cmd = f"samtools markdup -@ {threads} {sorted_bam} {final_bam}"
        run_command(markdup_cmd, "Samtools markdup (标记PCR重复)")
        
        # samtools markdup 不会自动创建索引，需要手动创建
        index_cmd = f"samtools index {final_bam}"
        run_command(index_cmd, "Samtools 索引")
    
    # 删除排序后的中间 BAM
    os.remove(sorted_bam)
    
    print(f"\n[完成] 样本 {sample_name} 处理完毕。输出: {final_bam}")
    return final_bam


def main():
    parser = argparse.ArgumentParser(
        description="将 FASTQ 数据比对到参考基因组并进行 BAM 处理（含 PCR 重复标记）。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    required = parser.add_argument_group("必需参数")
    required.add_argument(
        '-r', '--reference',
        required=True,
        metavar='REF.FASTA',
        help="参考基因组 FASTA 文件（需已建立索引）"
    )
    required.add_argument(
        '-p', '--plus-tag',
        required=True,
        metavar='SUFFIX',
        help="正向读段 (R1) 文件后缀，如: _1.fastq.gz"
    )
    required.add_argument(
        '-m', '--minus-tag',
        required=True,
        metavar='SUFFIX',
        help="反向读段 (R2) 文件后缀，如: _2.fastq.gz"
    )
    
    optional = parser.add_argument_group("可选参数")
    optional.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        metavar='N',
        help="使用的 CPU 线程数 (默认: 4)"
    )
    optional.add_argument(
        '-f', '--folder',
        default=".",
        metavar='DIR',
        help="FASTQ 文件所在文件夹 (默认: 当前目录)"
    )
    optional.add_argument(
        '-o', '--output',
        default=None,
        metavar='DIR',
        help="BAM 输出文件夹 (默认: 与输入文件夹相同)"
    )
    optional.add_argument(
        '--use-gatk-markdup',
        action='store_true',
        help="优先使用 GATK MarkDuplicates（更精确但更慢）"
    )
    
    args = parser.parse_args()
    
    # 处理路径
    input_folder = os.path.abspath(args.folder)
    output_folder = os.path.abspath(args.output) if args.output else input_folder
    reference = os.path.abspath(args.reference)
    
    # 检查参考基因组
    if not os.path.exists(reference):
        print(f"错误: 参考基因组不存在: {reference}")
        sys.exit(1)
    
    # 检查索引文件
    required_indices = [f"{reference}.fai", f"{reference}.bwt"]
    for idx_file in required_indices:
        if not os.path.exists(idx_file):
            print(f"错误: 缺少索引文件: {idx_file}")
            print("请先运行 1_index_reference.py 创建索引。")
            sys.exit(1)
    
    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)
    
    # 查找样本
    samples = find_samples(input_folder, args.plus_tag, args.minus_tag)
    
    if not samples:
        print(f"错误: 在 {input_folder} 中未找到匹配的 FASTQ 文件对。")
        print(f"  正向后缀: {args.plus_tag}")
        print(f"  反向后缀: {args.minus_tag}")
        sys.exit(1)
    
    print(f"\n找到 {len(samples)} 个样本:")
    for s in samples:
        print(f"  - {s['name']}")
    
    # 处理每个样本
    output_bams = []
    for sample in samples:
        bam_path = process_sample(
            sample, 
            reference, 
            args.threads, 
            output_folder,
            args.use_gatk_markdup
        )
        output_bams.append(bam_path)
    
    print("\n" + "="*60)
    print("所有样本处理完成！")
    print("="*60)
    print(f"输出 BAM 文件 ({len(output_bams)} 个):")
    for bam in output_bams:
        print(f"  - {bam}")
    print("\n下一步: 运行 3_gatk_call.py 进行变异检测。")


if __name__ == "__main__":
    main()