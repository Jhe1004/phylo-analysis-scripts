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
import logging
from datetime import datetime
from collections import OrderedDict


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


logger = setup_logger("6_combine_fasta")

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
        logger.error(f"读取文件 {os.path.basename(filepath)} 时出错: {e}")
        return None
    return sequences


def combine_fasta_files(output_filename="combined_sequences.fasta"):
    """
    将 consensus_fasta_output 目录中的所有 .fasta 文件合并成一个用于系统发育分析的 MSA 文件。
    """
    current_directory = os.getcwd()
    input_directory = os.path.join(current_directory, "consensus_fasta_output")
    
    if not os.path.exists(input_directory):
        logger.error(f"错误: 找不到输入目录 '{input_directory}'")
        return

    fasta_files_pattern = os.path.join(input_directory, "*.fasta")
    fasta_files = sorted(glob.glob(fasta_files_pattern)) # 排序以保证可重复性

    absolute_output_filename = os.path.join(current_directory, output_filename)

    if not fasta_files:
        logger.info(f"在目录 '{input_directory}' 中没有找到 .fasta 文件。")
        return

    logger.info(f"找到 {len(fasta_files)} 个 FASTA 文件待合并。")
    
    # --- 第一步：收集所有样本的 Scaffold 字典 ---
    all_sample_data = OrderedDict() # {sample_name: {scaffold_id: sequence}}
    scaffold_order = None # 用于记录第一个样本的 Scaffold 顺序作为标准

    for filepath in fasta_files:
        sample_name = os.path.splitext(os.path.basename(filepath))[0]
        if sample_name.endswith("_consensus"):
            sample_name = sample_name[:-len("_consensus")]
             
        logger.info(f"  正在解析: {os.path.basename(filepath)} (样本名: {sample_name})")
        
        parsed_seqs = parse_fasta(filepath)
        if parsed_seqs is None:
            logger.warning(f"  警告: 跳过无法解析的文件 {os.path.basename(filepath)}")
            continue
        
        if not parsed_seqs:
            logger.warning(f"  警告: 文件 {os.path.basename(filepath)} 中没有序列。跳过。")
            continue

        all_sample_data[sample_name] = parsed_seqs
        
        if scaffold_order is None:
            scaffold_order = list(parsed_seqs.keys())
            logger.info(f"  [INFO] 以此文件作为 Scaffold 顺序标准。检测到 {len(scaffold_order)} 个 Scaffolds/Contigs。")

    if not all_sample_data:
        logger.error("错误: 没有成功解析任何 FASTA 文件。")
        return
    
    if scaffold_order is None:
        logger.error("错误: 无法确定 Scaffold 顺序。")
        return

    # --- 第二步：检查所有样本的 Scaffold 是否一致 ---
    inconsistent_samples = []
    for sample_name, seqs in all_sample_data.items():
        sample_scaffolds = list(seqs.keys())
        if sample_scaffolds != scaffold_order:
            inconsistent_samples.append(sample_name)
            missing = set(scaffold_order) - set(sample_scaffolds)
            extra = set(sample_scaffolds) - set(scaffold_order)
            if missing:
                logger.warning(f"  警告: 样本 '{sample_name}' 缺少以下 Scaffolds: {missing}")
            if extra:
                logger.warning(f"  警告: 样本 '{sample_name}' 有多余的 Scaffolds: {extra}")

    if inconsistent_samples:
        logger.warning(f"\n警告: 以下 {len(inconsistent_samples)} 个样本的 Scaffold 组成与第一个样本不一致。")
        logger.warning("这可能导致合并后的矩阵不正确。建议检查原始数据或参考基因组。")

    # --- 第三步：合并并写入输出文件 ---
    logger.info(f"\n正在写入合并后的 FASTA 文件: {output_filename}")
    
    total_expected_length = None
    
    with open(absolute_output_filename, 'w') as outfile:
        for sample_name, seqs in all_sample_data.items():
            concatenated_seq_parts = []
            for scaffold_id in scaffold_order:
                if scaffold_id in seqs:
                    concatenated_seq_parts.append(seqs[scaffold_id])
                else:
                    logger.warning(f"  严重警告: 样本 '{sample_name}' 缺少 Scaffold '{scaffold_id}'，无法填充。")

            full_sequence = "".join(concatenated_seq_parts)
            
            if total_expected_length is None:
                total_expected_length = len(full_sequence)
            elif len(full_sequence) != total_expected_length:
                logger.warning(f"  警告: 样本 '{sample_name}' 的总序列长度 ({len(full_sequence)}) "
                      f"与预期 ({total_expected_length}) 不一致!")
             
            outfile.write(f">{sample_name}\n")
            line_width = 80
            for i in range(0, len(full_sequence), line_width):
                outfile.write(full_sequence[i:i+line_width] + "\n")

    logger.info(f"\n成功! 已将 {len(all_sample_data)} 个样本合并到 '{output_filename}'。")
    logger.info(f"每个样本的序列由 {len(scaffold_order)} 个 Scaffolds 拼接而成，总长度约 {total_expected_length:,} bp。")


if __name__ == "__main__":
    combine_fasta_files()