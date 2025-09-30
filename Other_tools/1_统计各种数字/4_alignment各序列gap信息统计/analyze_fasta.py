# -*- coding: utf-8 -*-

import os
import csv

# --- 请在这里修改参数 ---
# 请将 "your_alignment.fasta" 替换成您的FASTA文件名
fasta_file_path = "combined_sequences.fas"
# -------------------------


def analyze_fasta_to_csv(file_path):
    """
    解析FASTA文件，并将每条序列中非ATCG字符的比例分析结果保存到CSV文件中。

    Args:
        file_path (str): FASTA文件的路径。
    """
    # 检查文件是否存在
    if not os.path.exists(file_path):
        print(f"错误：文件 '{file_path}' 未找到。请检查文件名和路径是否正确。")
        return

    # --- 自动生成输出文件名 ---
    base_name = os.path.splitext(file_path)[0]
    csv_output_path = f"{base_name}_analysis.csv"
    
    # 定义标准碱基集合，用于快速查找（忽略大小写）
    standard_bases = {'A', 'T', 'C', 'G', 'a', 't', 'c', 'g'}
    
    sequences = {}
    current_header = None

    print(f"--- 开始分析文件: {file_path} ---")

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    current_header = line[1:]
                    sequences[current_header] = ""
                elif current_header:
                    sequences[current_header] += line
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return

    if not sequences:
        print("文件中没有找到任何FASTA格式的序列。")
        return

    # 准备用于写入CSV的数据
    results_data = []
    # CSV文件的表头
    csv_header = ["序列ID", "总长度(bp)", "非'ATCG'字符数", "所占比例(%)"]
    
    # 遍历所有解析出的序列并进行计算
    for header, sequence in sequences.items():
        total_length = len(sequence)
        
        if total_length == 0:
            non_standard_count = 0
            ratio_percent_str = "0.00"
        else:
            non_standard_count = sum(1 for base in sequence if base not in standard_bases)
            ratio = (non_standard_count / total_length) * 100
            ratio_percent_str = f"{ratio:.2f}" # 格式化为两位小数的字符串
        
        # 将结果添加到数据列表中
        results_data.append([header, total_length, non_standard_count, ratio_percent_str])

    # 将结果写入CSV文件
    try:
        with open(csv_output_path, 'w', newline='', encoding='utf-8-sig') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(csv_header)  # 写入表头
            writer.writerows(results_data) # 写入所有数据行
        
        print(f"\n分析完成！")
        print(f"结果已成功保存到文件: {csv_output_path}")

    except Exception as e:
        print(f"写入CSV文件时发生错误: {e}")


if __name__ == "__main__":
    # 直接调用主函数
    analyze_fasta_to_csv(fasta_file_path)