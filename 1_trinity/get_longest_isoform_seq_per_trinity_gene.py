#!/usr/bin/env python3

import sys
import re
import os
import glob
import argparse
from Bio import SeqIO

def process_fasta_file(input_filepath, output_filepath):
    """
    处理单个FASTA文件，找到每个基因的最长转录本并写入输出文件。

    Args:
        input_filepath (str): 输入的FASTA文件路径。
        output_filepath (str): 输出的FASTA文件路径。
    """
    # 用于存储每个基因最长转录本的字典
    gene_to_longest_transcript = {}
    
    # 编译用于从 Trinity ID 中提取基因 ID 的正则表达式
    gene_id_pattern = re.compile(r'^(.*(?:_g\d+|comp\d+_c\d+))_')

    # 读取并处理FASTA文件
    try:
        for record in SeqIO.parse(input_filepath, "fasta"):
            match = gene_id_pattern.search(record.id)
            
            if not match:
                print(
                    f"警告: 无法从序列ID '{record.id}' (文件: {os.path.basename(input_filepath)}) 中解析基因ID。跳过此记录。",
                    file=sys.stderr
                )
                continue

            gene_id = match.group(1)
            sequence_len = len(record.seq)

            # 检查是否需要更新该基因的最长转录本
            if (gene_id not in gene_to_longest_transcript or 
                    sequence_len > len(gene_to_longest_transcript[gene_id].seq)):
                gene_to_longest_transcript[gene_id] = record
    
    except FileNotFoundError:
        print(f"错误: 文件未找到 '{input_filepath}'", file=sys.stderr)
        return # 提前退出此文件的处理

    # 如果文件为空或未找到任何转录本，则直接返回
    if not gene_to_longest_transcript:
        print(f"信息: 在文件 '{os.path.basename(input_filepath)}' 中未找到有效转录本，不生成输出文件。", file=sys.stderr)
        return

    # 对找到的最长转录本按长度降序排序
    longest_transcripts = sorted(
        gene_to_longest_transcript.values(),
        key=lambda record: len(record.seq),
        reverse=True
    )

    # 将结果以 FASTA 格式写入对应的输出文件
    SeqIO.write(longest_transcripts, output_filepath, "fasta")
    print(f"-> 已处理 '{os.path.basename(input_filepath)}', 结果保存至 '{os.path.basename(output_filepath)}'")


def main():
    """
    主函数，负责解析参数和协调批量处理流程。
    """
    parser = argparse.ArgumentParser(
        description="批量处理一个文件夹内所有.fasta文件，为每个基因提取最长的转录本。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "input_directory",
        help="包含.fasta文件的输入文件夹路径"
    )
    args = parser.parse_args()

    # 打印警告信息
    print(
        "\n\n\tNOTE - longest transcript isn't always the best transcript! "
        "... consider filtering based on relative expression support ...\n\n",
        file=sys.stderr
    )
    
    input_dir = args.input_directory
    
    # 检查输入路径是否存在且为文件夹
    if not os.path.isdir(input_dir):
        print(f"错误: '{input_dir}' 不是一个有效的文件夹路径。", file=sys.stderr)
        sys.exit(1)

    # 使用 glob 查找所有 .fasta 后缀的文件
    # os.path.join 确保路径拼接的正确性
    fasta_files = glob.glob(os.path.join(input_dir, '*.fasta'))

    if not fasta_files:
        print(f"在文件夹 '{input_dir}' 中没有找到任何 .fasta 文件。", file=sys.stderr)
        sys.exit(0)

    print(f"共找到 {len(fasta_files)} 个 .fasta 文件进行处理...\n")

    # 遍历所有找到的fasta文件
    for infile_path in fasta_files:
        # 构建输出文件名：将 .fasta 替换为 .fas
        # os.path.splitext可以安全地分离文件名和扩展名
        base_name, _ = os.path.splitext(infile_path)
        outfile_path = base_name + '.fas'
        
        # 调用核心处理函数
        process_fasta_file(infile_path, outfile_path)
        
    print("\n\nok. 所有文件处理完毕.\n\n", file=sys.stderr)

if __name__ == "__main__":
    main()