#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from collections import defaultdict
from Bio import SeqIO

def parse_arguments():
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(description="统计指定物种在所有FASTA文件中出现的次数。")
    parser.add_argument(
        "-o", "--output_file",
        default="species_counts.txt",
        help="输出结果的文件名。默认是 species_counts.txt"
    )
    parser.add_argument(
        "-e", "--extensions",
        nargs="+",
        default=[".fasta", ".fa", ".fna", ".ffn", ".faa", ".fas"],
        help="要扫描的FASTA文件扩展名列表。默认包括 .fasta, .fa, .fna, .ffn, .faa, .frn"
    )
    return parser.parse_args()

def get_fasta_files(extensions):
    """
    获取当前目录下所有指定扩展名的FASTA文件
    """
    fasta_files = []
    for file in os.listdir("."):
        if any(file.lower().endswith(ext) for ext in extensions):
            fasta_files.append(file)
    return fasta_files

def extract_species_name(description):
    """
    从FASTA描述行中提取物种名称。
    假设物种名称是描述行的第一个部分，用空格或其他分隔符分开。
    根据实际情况调整此函数。
    例如：
    >SpeciesName1|GeneInfo -> SpeciesName1
    >SpeciesName2 GeneInfo -> SpeciesName2
    """
    # 使用空格或其他分隔符分割描述行
    if '|' in description:
        species = description.split('|')[0]
    else:
        species = description.split()[0]
    return species

def get_species_from_fasta(fasta_files):
    """
    自动从所有FASTA文件中提取物种名称，返回一个集合
    """
    species_set = set()
    for fasta in fasta_files:
        try:
            with open(fasta, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    description = record.description
                    species = extract_species_name(description)
                    species_set.add(species)
        except Exception as e:
            print(f"警告：无法读取文件 {fasta}。错误信息：{e}")
    return species_set

def count_species_in_fasta(fasta_files, species_set):
    """
    统计每个物种在FASTA文件中出现的次数
    """
    species_counts = defaultdict(int)
    for fasta in fasta_files:
        try:
            with open(fasta, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    description = record.description
                    species = extract_species_name(description)
                    if species in species_set:
                        species_counts[species] += 1
        except Exception as e:
            print(f"警告：无法读取文件 {fasta}。错误信息：{e}")
    return species_counts

def write_counts_to_file(species_counts, species_set, output_file):
    """
    将计数结果写入输出文件
    """
    try:
        with open(output_file, "w") as f:
            f.write("Species\tCount\n")
            for species in sorted(species_set):
                count = species_counts.get(species, 0)
                f.write(f"{species}\t{count}\n")
        print(f"计数结果已写入 {output_file}")
    except Exception as e:
        print(f"错误：无法写入文件 {output_file}。错误信息：{e}")
        sys.exit(1)

def main():
    args = parse_arguments()
    # 自动获取FASTA文件
    fasta_files = get_fasta_files(args.extensions)
    
    if not fasta_files:
        print("警告：未找到任何FASTA文件。请检查文件扩展名和当前目录。")
        sys.exit(1)
    
    # 从FASTA文件中提取物种名称
    species_set = get_species_from_fasta(fasta_files)
    
    # 统计物种出现次数
    species_counts = count_species_in_fasta(fasta_files, species_set)
    
    # 将计数结果写入文件
    write_counts_to_file(species_counts, species_set, args.output_file)

if __name__ == "__main__":
    main()