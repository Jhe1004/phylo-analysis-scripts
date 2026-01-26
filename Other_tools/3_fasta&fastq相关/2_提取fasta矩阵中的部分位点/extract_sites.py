#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import random
import sys

# --- 用户配置区 ---
# 请在这里修改您的参数

# 1. 输入文件：包含多序列比对的FASTA文件路径
INPUT_FASTA_FILE = "result.fasta"

# 2. 输出文件：提取部分位点后，新生成的FASTA文件路径
OUTPUT_FASTA_FILE = "extracted_sites.fasta"

# 3. 提取位点的百分比：输入一个0到100之间的数字
#    例如，输入 10.0 将随机提取10%的位点。
PERCENTAGE_TO_EXTRACT = 10.0

# --- 配置区结束 ---


def parse_fasta(filename: str) -> list[tuple[str, str]]:
    """
    高效且健壮地解析FASTA文件。
    可以处理序列中的换行符。
    返回一个包含 (header, sequence) 元组的列表，以保持原始顺序。
    (此函数与上一个脚本中的函数相同)
    """
    sequences = []
    current_seq_parts = []
    current_header = None

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header is not None:
                        sequences.append((current_header, "".join(current_seq_parts)))
                    
                    current_header = line[1:]
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
        
        if current_header is not None:
            sequences.append((current_header, "".join(current_seq_parts)))

    except FileNotFoundError:
        print(f"错误：输入文件未找到 -> {filename}", file=sys.stderr)
        sys.exit(1)

    if not sequences:
        print(f"错误：文件中没有找到任何FASTA格式的序列 -> {filename}", file=sys.stderr)
        sys.exit(1)
        
    # 为了兼容旧版Python，使用 'typing.List' 和 'typing.Tuple'
    # from typing import List, Tuple
    # -> List[Tuple[str, str]]
    return sequences

def main():
    """
    主执行函数
    """
    print("--- 开始运行FASTA位点提取脚本 ---")

    # 0. 验证配置
    if not 0.0 < PERCENTAGE_TO_EXTRACT <= 100.0:
        print(f"错误：提取百分比 ({PERCENTAGE_TO_EXTRACT}%) 必须在 (0, 100] 的范围内。", file=sys.stderr)
        sys.exit(1)
    
    # 1. 读取并解析输入的FASTA比对文件
    print(f"正在读取输入文件: {INPUT_FASTA_FILE}...")
    original_sequences = parse_fasta(INPUT_FASTA_FILE)
    
    # 2. 验证比对文件格式
    print("正在验证比对文件...")
    if not original_sequences:
        return 

    first_seq_len = len(original_sequences[0][1])
    for header, seq in original_sequences[1:]:
        if len(seq) != first_seq_len:
            print("错误：序列长度不一致！", file=sys.stderr)
            print(f"  序列 '{original_sequences[0][0]}' 的长度为 {first_seq_len}", file=sys.stderr)
            print(f"  但序列 '{header}' 的长度为 {len(seq)}", file=sys.stderr)
            print("请确保输入文件是合法的多序列比对文件。", file=sys.stderr)
            sys.exit(1)
            
    alignment_length = first_seq_len
    num_sequences = len(original_sequences)
    print(f"验证成功！共 {num_sequences} 条序列，原始比对长度为 {alignment_length}。")

    # 3. 计算需要抽取的位点数并进行随机抽样
    num_sites_to_extract = int(alignment_length * (PERCENTAGE_TO_EXTRACT / 100.0))
    
    if num_sites_to_extract == 0:
        print(f"警告：计算出的提取位点数为0（原始长度 {alignment_length}，百分比 {PERCENTAGE_TO_EXTRACT}%）。", file=sys.stderr)
        print("脚本已停止，不会生成空文件。", file=sys.stderr)
        sys.exit(0)

    print(f"将从 {alignment_length} 个位点中随机抽取 {num_sites_to_extract} 个位点 ({PERCENTAGE_TO_EXTRACT}%)。")

    # 核心步骤：从所有列的索引中，无放回地随机抽取指定数量的索引
    # random.sample 是实现此功能的完美函数
    all_indices = range(alignment_length)
    extracted_indices = sorted(random.sample(all_indices, k=num_sites_to_extract))

    # 4. 构建新的比对序列
    print("正在构建新的比对序列...")
    new_alignment = []
    for header, original_seq in original_sequences:
        # 使用列表推导式高效地从原始序列中挑出指定的列
        new_seq = "".join([original_seq[i] for i in extracted_indices])
        new_alignment.append((header, new_seq))

    # 5. 写入输出文件
    print(f"正在将结果写入文件: {OUTPUT_FASTA_FILE}...")
    try:
        # 确保输出目录存在
        output_dir = os.path.dirname(OUTPUT_FASTA_FILE)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            
        with open(OUTPUT_FASTA_FILE, 'w') as f_out:
            for header, new_seq in new_alignment:
                f_out.write(f">{header}\n")
                f_out.write(f"{new_seq}\n")
    except IOError as e:
        print(f"错误：无法写入文件 {OUTPUT_FASTA_FILE}。原因: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"新生成的比对文件长度为: {len(extracted_indices)}")
    print("--- 所有任务完成！ ---")


if __name__ == "__main__":
    # 确保在旧版Python上也能运行，我们注释掉了可能引发错误的类型提示
    # 对于Python 3.9+，可以保留 def parse_fasta(...) -> list[tuple[str, str]]:
    main()