# -*- coding: utf-8 -*-
"""
将多个样本的 Consensus FASTA 文件合并成一个多序列比对 (MSA) 矩阵。

本脚本专为处理包含多个 Scaffold/Contig 的参考基因组设计。
它会识别每个样本 FASTA 中的所有 Scaffold，并将来自不同样本的
相同 Scaffold 序列合并在一起（在同一个 Scaffold 条目下拼接）。

最终输出文件的格式为：
>Scaffold_1
SampleA_seq...SampleB_seq...SampleC_seq... (所有样本的 Scaffold_1 序列首尾相连)
>Scaffold_2
SampleA_seq...SampleB_seq...SampleC_seq...
...

或者，更常见和推荐的系统发育格式（每个样本一条总序列）：
>SampleA
Scaffold_1_seq...Scaffold_2_seq... (该样本所有 Scaffold 序列首尾相连)
>SampleB
Scaffold_1_seq...Scaffold_2_seq...
...

本脚本采用后者（按样本合并），更适合系统发育分析。
"""

import os
import glob
from collections import OrderedDict

def parse_fasta(filepath):
    """
    解析 FASTA 文件，返回一个 OrderedDict: {seq_id: sequence_string}
    使用 OrderedDict 以保留 Scaffold 的原始顺序。
    """
    sequences = OrderedDict()
    current_id = None
    current_seq_parts = []

    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    # 如果之前有序列，先保存
                    if current_id is not None:
                        sequences[current_id] = "".join(current_seq_parts)
                    # 处理新的序列头
                    # 序列 ID 一般是 '>' 之后到第一个空格之前的部分
                    current_id = line[1:].split()[0]
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
            # 保存最后一条序列
            if current_id is not None:
                sequences[current_id] = "".join(current_seq_parts)
    except Exception as e:
        print(f"读取文件 {os.path.basename(filepath)} 时出错: {e}")
        return None
    return sequences


def combine_fasta_files(output_filename="combined_sequences.fasta"):
    """
    将当前目录中的所有 .fasta 文件合并成一个用于系统发育分析的 MSA 文件。

    合并策略：
    - 每个样本的所有 Scaffold 序列会被**首尾拼接**成一条完整的序列。
    - 最终输出文件中的每条序列代表一个样本。
    - 拼接时，各样本内的 Scaffold 顺序会保持一致（按参考基因组中的顺序）。
    """
    current_directory = os.getcwd()
    fasta_files_pattern = os.path.join(current_directory, "*.fasta")
    fasta_files = sorted(glob.glob(fasta_files_pattern)) # 排序以保证可重复性

    absolute_output_filename = os.path.join(current_directory, output_filename)
    if absolute_output_filename in fasta_files:
        fasta_files.remove(absolute_output_filename)

    if not fasta_files:
        print("在当前目录中没有找到 .fasta 文件。")
        return

    print(f"找到 {len(fasta_files)} 个 FASTA 文件待合并。")
    
    # --- 第一步：收集所有样本的 Scaffold 字典 ---
    all_sample_data = OrderedDict() # {sample_name: {scaffold_id: sequence}}
    scaffold_order = None # 用于记录第一个样本的 Scaffold 顺序作为标准

    for filepath in fasta_files:
        sample_name = os.path.splitext(os.path.basename(filepath))[0]
        # 常见的后缀如 _consensus，可以去掉以简化样本名
        if sample_name.endswith("_consensus"):
            sample_name = sample_name[:-len("_consensus")]
            
        print(f"  正在解析: {os.path.basename(filepath)} (样本名: {sample_name})")
        
        parsed_seqs = parse_fasta(filepath)
        if parsed_seqs is None:
            print(f"  警告: 跳过无法解析的文件 {os.path.basename(filepath)}")
            continue
        
        if not parsed_seqs:
            print(f"  警告: 文件 {os.path.basename(filepath)} 中没有序列。跳过。")
            continue

        all_sample_data[sample_name] = parsed_seqs
        
        # 记录第一个成功解析的样本的 Scaffold 顺序
        if scaffold_order is None:
            scaffold_order = list(parsed_seqs.keys())
            print(f"  [INFO] 以此文件作为 Scaffold 顺序标准。检测到 {len(scaffold_order)} 个 Scaffolds/Contigs。")

    if not all_sample_data:
        print("错误: 没有成功解析任何 FASTA 文件。")
        return
    
    if scaffold_order is None:
        print("错误: 无法确定 Scaffold 顺序。")
        return

    # --- 第二步：检查所有样本的 Scaffold 是否一致 ---
    inconsistent_samples = []
    for sample_name, seqs in all_sample_data.items():
        sample_scaffolds = list(seqs.keys())
        if sample_scaffolds != scaffold_order:
            inconsistent_samples.append(sample_name)
            # 进行更详细的检查
            missing = set(scaffold_order) - set(sample_scaffolds)
            extra = set(sample_scaffolds) - set(scaffold_order)
            if missing:
                print(f"  警告: 样本 '{sample_name}' 缺少以下 Scaffolds: {missing}")
            if extra:
                print(f"  警告: 样本 '{sample_name}' 有多余的 Scaffolds: {extra}")

    if inconsistent_samples:
        print(f"\n警告: 以下 {len(inconsistent_samples)} 个样本的 Scaffold 组成与第一个样本不一致。")
        print("这可能导致合并后的矩阵不正确。建议检查原始数据或参考基因组。")
        # 可以选择是否继续，这里选择继续但给出警告

    # --- 第三步：合并并写入输出文件 ---
    # 合并策略：每个样本，将其所有 Scaffold 按照 scaffold_order 拼接成一条长序列
    print(f"\n正在写入合并后的 FASTA 文件: {output_filename}")
    
    total_expected_length = None
    
    with open(absolute_output_filename, 'w') as outfile:
        for sample_name, seqs in all_sample_data.items():
            concatenated_seq_parts = []
            for scaffold_id in scaffold_order:
                if scaffold_id in seqs:
                    concatenated_seq_parts.append(seqs[scaffold_id])
                else:
                    # 如果某个样本缺少某个 Scaffold，用 N 填充（需要知道长度）
                    # 这里简化处理，假设不会缺失。真实情况应从其他样本获取长度。
                    print(f"  严重警告: 样本 '{sample_name}' 缺少 Scaffold '{scaffold_id}'，无法填充。")
                    # 此处可根据需要添加填充逻辑

            full_sequence = "".join(concatenated_seq_parts)
            
            if total_expected_length is None:
                total_expected_length = len(full_sequence)
            elif len(full_sequence) != total_expected_length:
                print(f"  警告: 样本 '{sample_name}' 的总序列长度 ({len(full_sequence)}) "
                      f"与预期 ({total_expected_length}) 不一致!")
            
            outfile.write(f">{sample_name}\n")
            # 写入序列，每行 80 个字符以提高可读性
            line_width = 80
            for i in range(0, len(full_sequence), line_width):
                outfile.write(full_sequence[i:i+line_width] + "\n")

    print(f"\n成功! 已将 {len(all_sample_data)} 个样本合并到 '{output_filename}'。")
    print(f"每个样本的序列由 {len(scaffold_order)} 个 Scaffolds 拼接而成，总长度约 {total_expected_length:,} bp。")


if __name__ == "__main__":
    combine_fasta_files()