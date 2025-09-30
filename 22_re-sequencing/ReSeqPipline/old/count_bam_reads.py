#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import pysam
from collections import defaultdict

def parse_arguments():
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(description="统计当前文件夹中所有 BAM 文件的 reads 数量。")
    parser.add_argument(
        "-o", "--output",
        default="bam_read_counts.txt",
        help="输出结果的文件名。默认为 bam_read_counts.txt"
    )
    parser.add_argument(
        "-e", "--extensions",
        nargs="+",
        default=[".bam"],
        help="要扫描的 BAM 文件扩展名列表。默认是 .bam"
    )
    return parser.parse_args()

def get_bam_files(extensions):
    """
    获取当前目录下所有指定扩展名的 BAM 文件
    """
    bam_files = []
    for file in os.listdir("."):
        if any(file.lower().endswith(ext) for ext in extensions):
            bam_files.append(file)
    return bam_files

def check_index_file(bam_file):
    """
    检查 BAM 文件是否有对应的 .csi 索引文件
    """
    index_file = bam_file + ".csi"
    if os.path.isfile(index_file):
        return True
    else:
        print(f"警告：索引文件 {index_file} 未找到。跳过 {bam_file}。")
        return False

def count_reads(bam_file):
    """
    使用 pysam 统计 BAM 文件中的 reads 数量
    """
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            read_count = bam.count()
        return read_count
    except Exception as e:
        print(f"错误：无法读取 BAM 文件 {bam_file}。错误信息：{e}")
        return None

def write_counts_to_file(counts_dict, output_file):
    """
    将统计结果写入输出文件
    """
    try:
        with open(output_file, "w") as f:
            f.write("BAM_File\tRead_Count\n")
            for bam, count in counts_dict.items():
                f.write(f"{bam}\t{count}\n")
        print(f"\n计数结果已写入 {output_file}")
    except Exception as e:
        print(f"错误：无法写入文件 {output_file}。错误信息：{e}")
        sys.exit(1)

def main():
    args = parse_arguments()
    bam_files = get_bam_files(args.extensions)

    if not bam_files:
        print("警告：未找到任何 BAM 文件。请检查文件扩展名和当前目录。")
        sys.exit(0)

    counts_dict = defaultdict(int)

    print(f"找到 {len(bam_files)} 个 BAM 文件。开始统计 reads 数量...\n")
    
    for bam in bam_files:
        if check_index_file(bam):
            read_count = count_reads(bam)
            if read_count is not None:
                counts_dict[bam] = read_count
                print(f"{bam}: {read_count} reads")
        else:
            print(f"跳过 {bam}。\n")

    if counts_dict:
        write_counts_to_file(counts_dict, args.output)
    else:
        print("没有有效的 BAM 文件被处理。")

if __name__ == "__main__":
    main()