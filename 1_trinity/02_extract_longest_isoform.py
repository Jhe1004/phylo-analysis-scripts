# -*- coding: utf-8 -*-
"""
脚本名称: 02_extract_longest_isoform.py
功能: 遍历文件夹内的 FASTA 文件，提取每个基因的最长转录本
作者: Dr. Jian He (Refactored by Gemini)
依赖: Python 3.6+, Biopython (pip install biopython)
"""

import os
import sys
import re
import glob

# 尝试导入 Biopython，如果没有安装则给出提示
try:
    from Bio import SeqIO
except ImportError:
    print("错误: 缺少 'Biopython' 库。")
    print("请运行: pip install biopython")
    sys.exit(1)

# ======================= ### CONFIGURATION (用户配置区) ### =======================
# 1. 路径设置
# 包含 FASTA 文件的输入文件夹路径
INPUT_DIRECTORY = "."  

# 2. 文件扩展名设置
# 需要处理的输入文件后缀 (例如 "*.fasta" 或 "*.Trinity.fasta")
INPUT_FILE_PATTERN = "*.fasta"

# 输出文件的后缀名 (将在原文件名基础上去除后缀并添加此后缀)
OUTPUT_EXTENSION = ".longest.fas"

# 3. 基因ID解析设置
# 用于从序列 ID 中提取基因 ID 的正则表达式
# Trinity 默认格式示例: TRINITY_DN1000_c0_g1_i1 -> 基因ID: TRINITY_DN1000_c0_g1
# 默认正则捕获最后一个下划线之前的所有内容
GENE_ID_REGEX = r'^(.*)_i\d+$' 
# =================================================================================


def get_gene_id(seq_id, pattern):
    """
    使用正则表达式从序列 ID 中解析基因 ID。
    """
    match = pattern.search(seq_id)
    if match:
        return match.group(1)
    return None

def process_single_file(filepath, output_path, regex_pattern):
    """
    读取单个 FASTA 文件，筛选最长转录本并写入新文件。
    """
    gene_dict = {} # 格式: { 'gene_id': SeqRecord_object }
    
    print(f"正在读取: {os.path.basename(filepath)} ...")
    
    count_total = 0
    count_kept = 0
    
    try:
        # 第一遍扫描: 找出每个基因的最长序列
        for record in SeqIO.parse(filepath, "fasta"):
            count_total += 1
            gene_id = get_gene_id(record.id, regex_pattern)
            
            if not gene_id:
                # 如果正则匹配失败，尝试使用稍微宽松的 Trinity 默认逻辑 (兼容旧版脚本逻辑)
                # 旧逻辑: split('_') 然后去掉最后一部分
                parts = record.id.rsplit('_', 1)
                if len(parts) > 1:
                    gene_id = parts[0]
                else:
                    print(f"  [警告] 无法解析 ID: {record.id}，跳过。")
                    continue
            
            current_len = len(record.seq)
            
            # 核心比较逻辑
            if gene_id not in gene_dict:
                gene_dict[gene_id] = record
            else:
                existing_len = len(gene_dict[gene_id].seq)
                if current_len > existing_len:
                    gene_dict[gene_id] = record
        
        # 第二遍: 写入结果 (按长度降序排列，美观且实用)
        longest_transcripts = sorted(
            gene_dict.values(), 
            key=lambda x: len(x.seq), 
            reverse=True
        )
        
        if longest_transcripts:
            count_kept = len(longest_transcripts)
            SeqIO.write(longest_transcripts, output_path, "fasta")
            print(f"  -> 完成。输入序列: {count_total}, 提取基因: {count_kept}")
            print(f"  -> 结果已保存至: {output_path}\n")
        else:
            print("  [警告] 该文件中未找到有效序列，不生成输出文件。\n")
            
    except Exception as e:
        print(f"  [错误] 处理文件 {filepath} 时发生异常: {e}\n")

def main():
    # 检查输入目录
    if not os.path.isdir(INPUT_DIRECTORY):
        print(f"错误: 输入目录不存在 - {INPUT_DIRECTORY}")
        sys.exit(1)

    # 编译正则对象
    gene_pattern = re.compile(GENE_ID_REGEX)
    
    # 查找文件
    search_path = os.path.join(INPUT_DIRECTORY, INPUT_FILE_PATTERN)
    files = glob.glob(search_path)
    
    if not files:
        print(f"在目录 '{INPUT_DIRECTORY}' 中未找到匹配 '{INPUT_FILE_PATTERN}' 的文件。")
        sys.exit(0)
        
    print(f"找到 {len(files)} 个文件，开始处理...\n")
    print("-" * 30)

    for infile in files:
        # 构建输出文件名
        # os.path.splitext 只会去掉最后一个扩展名，比如 'sample.fasta' -> 'sample'
        file_base = os.path.splitext(infile)[0]
        outfile = file_base + OUTPUT_EXTENSION
        
        process_single_file(infile, outfile, gene_pattern)
        
    print("-" * 30)
    print("所有任务处理完毕。")

if __name__ == "__main__":
    main()