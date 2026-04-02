#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
为参考基因组创建所有必需的索引文件。

本脚本会依次执行以下命令：
1. bwa index: 创建 BWA 比对索引 (.amb, .ann, .bwt, .pac, .sa)
2. samtools faidx: 创建 FASTA 索引 (.fai)
3. gatk CreateSequenceDictionary: 创建序列字典 (.dict) - GATK 必须

软件依赖:
- bwa
- samtools
- gatk (或 picard)
"""

import os
import argparse
import subprocess
import sys
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


logger = setup_logger("1_index_reference")


def run_command(command, description):
    """运行命令并记录日志"""
    logger.info(f"\n{'='*60}")
    logger.info(f"[步骤] {description}")
    logger.debug(f"[命令] {command}")
    logger.info('='*60)
    
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        logger.error(f"[错误] 命令执行失败，返回码: {result.returncode}")
        sys.exit(1)
    logger.info(f"[完成] {description}")


def index_reference(reference_path):
    """为参考基因组创建所有索引"""
    
    # 检查参考基因组是否存在
    if not os.path.exists(reference_path):
        logger.error(f"参考基因组文件不存在: {reference_path}")
        sys.exit(1)
    
    ref_dir = os.path.dirname(os.path.abspath(reference_path))
    ref_basename = os.path.basename(reference_path)
    ref_name_no_ext = os.path.splitext(ref_basename)[0]
    
    logger.info(f"\n参考基因组: {reference_path}")
    logger.info(f"索引文件将保存在: {ref_dir}")
    
    # 1. BWA 索引
    # 对于大于 2GB 的基因组使用 bwtsw 算法，否则使用 is 算法
    # 这里默认使用 bwtsw，适用于所有大小
    bwa_cmd = f"bwa index -a bwtsw {reference_path}"
    run_command(bwa_cmd, "创建 BWA 索引")
    
    # 2. Samtools FASTA 索引
    faidx_cmd = f"samtools faidx {reference_path}"
    run_command(faidx_cmd, "创建 Samtools FASTA 索引 (.fai)")
    
    # 3. GATK 序列字典 (.dict)
    dict_path = os.path.join(ref_dir, f"{ref_name_no_ext}.dict")
    # 如果已存在则先删除（GATK 不会覆盖）
    if os.path.exists(dict_path):
        logger.info(f"[信息] 发现已存在的 .dict 文件，将删除后重新生成: {dict_path}")
        os.remove(dict_path)
    
    # 尝试使用 gatk，如果失败则尝试 picard
    gatk_cmd = f"gatk CreateSequenceDictionary -R {reference_path} -O {dict_path}"
    try:
        run_command(gatk_cmd, "创建 GATK 序列字典 (.dict)")
    except SystemExit:
        logger.info("[信息] GATK 命令失败，尝试使用 Picard...")
        picard_cmd = f"picard CreateSequenceDictionary R={reference_path} O={dict_path}"
        run_command(picard_cmd, "创建 Picard 序列字典 (.dict)")
    
    logger.info("\n" + "="*60)
    logger.info("所有索引文件创建完成！")
    logger.info("="*60)
    logger.info(f"  - BWA 索引: {reference_path}.* (多个文件)")
    logger.info(f"  - FASTA 索引: {reference_path}.fai")
    logger.info(f"  - 序列字典: {dict_path}")


def main():
    parser = argparse.ArgumentParser(
        description="为参考基因组创建 BWA, Samtools 和 GATK 所需的所有索引文件。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-r', '--reference',
        required=True,
        metavar='REF.FASTA',
        help="参考基因组 FASTA 文件路径"
    )
    
    args = parser.parse_args()
    index_reference(args.reference)


if __name__ == "__main__":
    main()