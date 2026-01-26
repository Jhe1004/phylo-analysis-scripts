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


def run_command(command, description):
    """运行命令并打印状态"""
    print(f"\n{'='*60}")
    print(f"[步骤] {description}")
    print(f"[命令] {command}")
    print('='*60)
    
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print(f"[错误] 命令执行失败，返回码: {result.returncode}")
        sys.exit(1)
    print(f"[完成] {description}")


def index_reference(reference_path):
    """为参考基因组创建所有索引"""
    
    # 检查参考基因组是否存在
    if not os.path.exists(reference_path):
        print(f"错误: 参考基因组文件不存在: {reference_path}")
        sys.exit(1)
    
    ref_dir = os.path.dirname(os.path.abspath(reference_path))
    ref_basename = os.path.basename(reference_path)
    ref_name_no_ext = os.path.splitext(ref_basename)[0]
    
    print(f"\n参考基因组: {reference_path}")
    print(f"索引文件将保存在: {ref_dir}")
    
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
        print(f"[信息] 发现已存在的 .dict 文件，将删除后重新生成: {dict_path}")
        os.remove(dict_path)
    
    # 尝试使用 gatk，如果失败则尝试 picard
    gatk_cmd = f"gatk CreateSequenceDictionary -R {reference_path} -O {dict_path}"
    try:
        run_command(gatk_cmd, "创建 GATK 序列字典 (.dict)")
    except SystemExit:
        print("[信息] GATK 命令失败，尝试使用 Picard...")
        picard_cmd = f"picard CreateSequenceDictionary R={reference_path} O={dict_path}"
        run_command(picard_cmd, "创建 Picard 序列字典 (.dict)")
    
    print("\n" + "="*60)
    print("所有索引文件创建完成！")
    print("="*60)
    print(f"  - BWA 索引: {reference_path}.* (多个文件)")
    print(f"  - FASTA 索引: {reference_path}.fai")
    print(f"  - 序列字典: {dict_path}")


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