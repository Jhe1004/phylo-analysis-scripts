#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import random
import sys
from collections import defaultdict
from typing import List, Tuple  # <--- 第1步：在这里添加导入

# --- 用户配置区 ---
# 请在这里修改您的参数

# 1. 输入文件：包含多序列比对的FASTA文件路径
INPUT_FASTA_FILE = "result.fasta"

# 2. 输出目录：用于存放生成的bootstrap文件的文件夹
OUTPUT_DIR = "bootstrap_output"

# 3. 输出文件前缀：生成的文件将命名为 "前缀_1.fasta", "前缀_2.fasta", ...
OUTPUT_PREFIX = "bootstrap_replicate"

# 4. Bootstrap重抽样次数：您希望生成多少个新的比对文件
NUM_REPLICATES = 10

# --- 配置区结束 ---


# <--- 第2步：修改下面这一行
def parse_fasta(filename: str) -> List[Tuple[str, str]]:
    """
    高效且健壮地解析FASTA文件。
    可以处理序列中的换行符。
    返回一个包含 (header, sequence) 元组的列表，以保持原始顺序。
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
                    # 如果已有记录，先保存上一个序列
                    if current_header is not None:
                        sequences.append((current_header, "".join(current_seq_parts)))
                    
                    # 记录新的序列头，并重置序列部分
                    current_header = line[1:]
                    current_seq_parts = []
                else:
                    # 添加序列部分
                    current_seq_parts.append(line)
        
        # 添加文件中的最后一个序列
        if current_header is not None:
            sequences.append((current_header, "".join(current_seq_parts)))

    except FileNotFoundError:
        print(f"错误：输入文件未找到 -> {filename}", file=sys.stderr)
        sys.exit(1) # 脚本异常退出

    if not sequences:
        print(f"错误：文件中没有找到任何FASTA格式的序列 -> {filename}", file=sys.stderr)
        sys.exit(1)
        
    return sequences

def main():
    """
    主执行函数
    """
    print("--- 开始运行Bootstrap脚本 ---")
    
    # 1. 读取并解析输入的FASTA比对文件
    print(f"正在读取输入文件: {INPUT_FASTA_FILE}...")
    original_sequences = parse_fasta(INPUT_FASTA_FILE)
    
    # 2. 健壮性检查：验证所有序列长度是否一致
    print("正在验证比对文件...")
    if not original_sequences:
        # 这个检查在 `parse_fasta` 中已经包含了，但为了清晰再加一次
        return 

    first_seq_len = len(original_sequences[0][1])
    for header, seq in original_sequences[1:]:
        if len(seq) != first_seq_len:
            print(f"错误：序列长度不一致！", file=sys.stderr)
            print(f"  序列 '{original_sequences[0][0]}' 的长度为 {first_seq_len}", file=sys.stderr)
            print(f"  但序列 '{header}' 的长度为 {len(seq)}", file=sys.stderr)
            print("请确保输入文件是合法的多序列比对文件。", file=sys.stderr)
            sys.exit(1)
            
    alignment_length = first_seq_len
    num_sequences = len(original_sequences)
    print(f"验证成功！共 {num_sequences} 条序列，比对长度为 {alignment_length}。")

    # 3. 准备输出目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"输出文件将保存在目录: {OUTPUT_DIR}")

    # 4. 执行Bootstrap重抽样
    print(f"开始生成 {NUM_REPLICATES} 个Bootstrap重抽样文件...")
    
    # 生成一个从0到alignment_length-1的索引列表，用于抽样
    column_indices = list(range(alignment_length))

    for i in range(1, NUM_REPLICATES + 1):
        # 核心步骤：有放回地随机抽取列索引
        # random.choices 是最高效的实现方式
        bootstrap_indices = random.choices(column_indices, k=alignment_length)

        # 构建新的比对文件
        output_filename = os.path.join(OUTPUT_DIR, f"{OUTPUT_PREFIX}_{i}.fasta")
        
        try:
            with open(output_filename, 'w') as f_out:
                # 对原始比对中的每一条序列，根据新的列索引构建新序列
                for header, original_seq in original_sequences:
                    # 使用列表推导式和join方法，这是Python中构建字符串最快的方式
                    new_seq = "".join([original_seq[j] for j in bootstrap_indices])
                    
                    # 写入输出文件，格式为 >header\nsequence\n
                    f_out.write(f">{header}\n")
                    f_out.write(f"{new_seq}\n")
        except IOError as e:
            print(f"错误：无法写入文件 {output_filename}。原因: {e}", file=sys.stderr)
            sys.exit(1)

        print(f"  已生成文件: {output_filename}")

    print("--- 所有任务完成！ ---")


if __name__ == "__main__":
    # 当脚本被直接执行时，运行main函数
    main()