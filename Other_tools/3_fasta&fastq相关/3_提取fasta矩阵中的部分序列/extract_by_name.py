#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

# --- 用户配置区 ---
# 请在这里修改您的参数

# 1. 输入的主FASTA文件：包含所有序列的矩阵文件
FASTA_INPUT_FILE = "all_sequences.fasta"

# 2. 名称列表文件：一个文本文件，每行包含一个需要提取的序列名
NAMES_FILE = "names_to_extract.txt"

# 3. 输出的FASTA文件：存放提取出的序列
FASTA_OUTPUT_FILE = "extracted_sequences.fasta"

# --- 配置区结束 ---


def read_names_to_extract(filename: str) -> set:
    """
    从一个文本文件中读取序列名称。
    返回一个包含所有名称的集合(set)，以便快速查找。
    """
    names = set()
    try:
        with open(filename, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                if stripped_line:  # 忽略空行
                    names.add(stripped_line)
    except FileNotFoundError:
        print(f"错误：名称列表文件未找到 -> {filename}", file=sys.stderr)
        sys.exit(1)
    
    if not names:
        print(f"警告：名称列表文件 '{filename}' 为空或未能读取到任何名称。", file=sys.stderr)
    
    return names


def parse_fasta_generator(filename: str):
    """
    一个内存高效的FASTA解析器，使用生成器(generator)逐条产生序列。
    这样可以处理任意大小的FASTA文件而不会耗尽内存。
    产出 (header, sequence) 元组。
    """
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
                        yield current_header, "".join(current_seq_parts)
                    
                    current_header = line[1:]
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
        
        if current_header is not None:
            yield current_header, "".join(current_seq_parts)
            
    except FileNotFoundError:
        print(f"错误：主FASTA文件未找到 -> {filename}", file=sys.stderr)
        sys.exit(1)


def main():
    """
    主执行函数
    """
    print("--- 开始根据名称提取FASTA序列 ---")

    # 1. 读取需要提取的序列名称
    print(f"正在读取名称列表: {NAMES_FILE}...")
    names_to_extract = read_names_to_extract(NAMES_FILE)
    total_names_requested = len(names_to_extract)
    if total_names_requested == 0:
        print("由于未提供任何需要提取的名称，脚本已停止。")
        sys.exit(0)
    print(f"读取完成，共需要提取 {total_names_requested} 个序列。")

    # 2. 准备一个集合用于追踪未找到的名称
    missing_names = names_to_extract.copy()
    found_count = 0

    print(f"正在扫描主FASTA文件: {FASTA_INPUT_FILE}...")
    
    try:
        with open(FASTA_OUTPUT_FILE, 'w') as f_out:
            # 3. 逐条解析FASTA文件并进行匹配
            for header, sequence in parse_fasta_generator(FASTA_INPUT_FILE):
                # 智能匹配：通常只比较'>'后第一个空格前的内容
                # 例如，匹配 'gene1' 到 '>gene1 description...'
                identifier = header.split()[0]

                if identifier in names_to_extract:
                    f_out.write(f">{header}\n")
                    f_out.write(f"{sequence}\n")
                    found_count += 1
                    # 从待查找集合中移除已找到的名称
                    missing_names.discard(identifier)
                    
                    # 优化：如果所有需要的序列都找到了，可以提前结束
                    if not missing_names:
                        break
    
    except IOError as e:
        print(f"错误：无法写入输出文件 {FASTA_OUTPUT_FILE}。原因: {e}", file=sys.stderr)
        sys.exit(1)

    # 4. 生成总结报告
    print("\n--- 任务总结 ---")
    print(f"成功写入 {found_count} 条序列到文件: {FASTA_OUTPUT_FILE}")

    if missing_names:
        print(f"\n警告：有 {len(missing_names)} 个请求的序列名未在输入文件中找到：")
        for name in sorted(list(missing_names)):
            print(f"  - {name}")
    else:
        print("所有请求的序列都已成功找到并提取！")
    
    print("\n--- 所有任务完成！ ---")


if __name__ == "__main__":
    main()