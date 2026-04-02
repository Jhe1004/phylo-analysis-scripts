#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 FASTQ 测序数据比对到参考基因组，并进行标准的 BAM 处理。

本脚本执行以下步骤：
1. BWA MEM: 将双端测序数据比对到参考基因组
2. Samtools sort: 对 BAM 文件按坐标排序
3. Samtools index: 为最终 BAM 文件建立索引

软件依赖:
- bwa
- samtools
"""

import os
import argparse
import subprocess
import sys
import glob
import concurrent.futures
import logging
from datetime import datetime


def setup_logger(script_name):
    """配置日志记录器，同时输出到文件和控制台"""
    log_file = f"{script_name}.log"
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.DEBUG)
    
    logger.handlers.clear()
    
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    logger.info(f"日志文件: {log_file}")
    return logger


logger = setup_logger("2_bwa_map")


def run_command(command, description, sample_name, exit_on_fail=True):
    """运行命令并记录日志，带有样本名标识"""
    prefix = f"[{sample_name}]"
    logger.info(f"\n{prefix} [子步骤] {description}")
    logger.debug(f"{prefix} [命令] {command}")
    
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        logger.error(f"{prefix} [错误] 命令执行失败，返回码: {result.returncode}")
        if exit_on_fail:
            raise RuntimeError(f"Sample {sample_name} failed at step: {description}")
        return False
    return True


def find_samples(folder, plus_tag, minus_tag):
    """在指定文件夹中查找所有样本，支持多种常见后缀格式"""
    samples_dict = {}
    
    # 定义所有可能的后缀对
    tag_pairs = [
        ("_1.clean.fq.gz", "_2.clean.fq.gz"),
        ("_1.fq.gz", "_2.fq.gz"),
        ("_R1.fastq.gz", "_R2.fastq.gz"),
        ("_R1.fq.gz", "_R2.fq.gz")
    ]
    
    # 优先使用用户提供的 tag，将其放在列表首位
    if (plus_tag, minus_tag) not in tag_pairs:
        tag_pairs.insert(0, (plus_tag, minus_tag))

    for p_tag, m_tag in tag_pairs:
        for filepath in glob.glob(os.path.join(folder, f"*{p_tag}")):
            filename = os.path.basename(filepath)
            sample_name = filename.replace(p_tag, "")
            
            # 如果该样本已经通过之前的 tag 找到了，则跳过
            if sample_name in samples_dict:
                continue
                
            minus_file = os.path.join(folder, f"{sample_name}{m_tag}")
            if os.path.exists(minus_file):
                samples_dict[sample_name] = {
                    'name': sample_name,
                    'r1': filepath,
                    'r2': minus_file
                }
    
    if not samples_dict:
        return []
        
    logger.info(f"匹配到 {len(samples_dict)} 个样本。")
    return list(samples_dict.values())


def process_sample(sample, reference, threads, output_folder):
    """处理单个样本的完整比对流程"""
    
    sample_name = sample['name']
    r1_path = sample['r1']
    r2_path = sample['r2']
    
    # 定义进度输出前缀
    p = f"[{sample_name}]"
    
    logger.info(f"\n{p} {'='*40}")
    logger.info(f"{p} 开始处理: R1={os.path.basename(r1_path)}, R2={os.path.basename(r2_path)}")
    logger.info(f"{p} {'='*40}")
    
    # 定义中间文件和输出文件路径
    raw_bam = os.path.join(output_folder, f"{sample_name}.raw.bam")
    sorted_bam = os.path.join(output_folder, f"{sample_name}.sorted.bam")
    final_bam = os.path.join(output_folder, f"{sample_name}.bam")
    
    # 断点续传检查：如果最终 BAM 已存在，直接跳过
    if os.path.exists(final_bam) and os.path.getsize(final_bam) > 0:
        logger.info(f"{p} 最终 BAM 已存在，跳过。")
        return final_bam
    
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
        
        # 步骤 3: 直接重命名排序后的BAM为最终BAM (跳过MarkDuplicates)
        if os.path.exists(sorted_bam):
            os.rename(sorted_bam, final_bam)
        
        # 步骤 4: 索引最终BAM
        index_cmd = f"samtools index {final_bam}"
        run_command(index_cmd, "Samtools 索引", sample_name)
        
        logger.info(f"\n{p} [完成] 处理完毕。输出: {final_bam}")
        return final_bam
    except Exception as e:
        logger.error(f"\n{p} [严重错误] 样本处理中断: {e}")
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
    
    args = parser.parse_args()
    
    input_folder = os.path.abspath(args.folder)
    output_folder = os.path.abspath(args.output) if args.output else input_folder
    reference = os.path.abspath(args.reference)
    
    # 检查环境和文件... (省略部分检查)
    os.makedirs(output_folder, exist_ok=True)
    samples = find_samples(input_folder, args.plus_tag, args.minus_tag)
    
    if not samples:
        logger.error("未找到有效样本，退出。")
        sys.exit(1)

    num_parallel = min(args.parallel_samples, len(samples))
    threads_per_sample = args.threads
    total_threads = num_parallel * threads_per_sample
    
    logger.info(f"\n[任务启动] 样本总数: {len(samples)}")
    logger.info(f"[资源配置] 并行样本数: {num_parallel}")
    logger.info(f"[资源配置] 每样本线程: {threads_per_sample}")
    logger.info(f"[资源配置] 预计峰值总线程: {total_threads}")
    
    output_bams = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel) as executor:
        future_to_sample = {
            executor.submit(
                process_sample, 
                sample, 
                reference, 
                threads_per_sample, 
                output_folder
            ): sample['name'] for sample in samples
        }
        
        for future in concurrent.futures.as_completed(future_to_sample):
            sample_name = future_to_sample[future]
            try:
                bam_path = future.result()
                if bam_path:
                    output_bams.append(bam_path)
            except Exception as exc:
                logger.error(f"{sample_name} 生成了未预料的异常: {exc}")

    logger.info("\n" + "="*60)
    logger.info(f"所有任务处理尝试完成！成功: {len(output_bams)}/{len(samples)}")
    logger.info("="*60)
    for bam in output_bams:
        logger.info(f"  - {bam}")


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()