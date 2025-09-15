# -*- coding: utf-8 -*-

import csv
import sys

# 尝试导入Biopython库，如果失败则给出提示
try:
    from Bio import SeqIO
except ImportError:
    print("错误: 本脚本需要 Biopython 库。")
    print("请使用 'pip install biopython' 命令进行安装。")
    sys.exit(1)

# =============================================================================
#                            用户配置区域
#
#  请在这里修改你的输入和输出文件路径。
#  你只需要更改引号内的文件名即可。
# =============================================================================

# 1. 输入的FASTA文件路径 (包含原始序列的文件)
FASTA_FILE_PATH = 'input.fasta'

# 2. 名称映射CSV文件路径
#    CSV文件格式:
#    第1列: 原始序列名 (FASTA文件中'>'后的ID)
#    第2列: 你想要替换成的新序列名
#    注意: CSV文件不应包含表头(header)。
MAPPING_FILE_PATH = 'mapping.csv'

# 3. 输出的FASTA文件路径 (结果将保存到这里)
OUTPUT_FASTA_PATH = 'output_renamed_linear.fasta'


# =============================================================================
#                            脚本核心逻辑
#
#            一般来说，你不需要修改下面的代码。
# =============================================================================

def create_replacement_map(mapping_file):
    """
    从CSV文件中读取名称映射关系，并创建一个字典。
    """
    replacement_dict = {}
    try:
        with open(mapping_file, mode='r', encoding='utf-8') as infile:
            reader = csv.reader(infile)
            for rows in reader:
                if len(rows) >= 2:
                    old_name = rows[0].strip()
                    new_name = rows[1].strip()
                    replacement_dict[old_name] = new_name
                else:
                    print(f"警告: 跳过格式不正确的行: {rows}")
        print(f"成功: 从 '{mapping_file}' 文件中加载了 {len(replacement_dict)} 个名称映射关系。")
        return replacement_dict
    except FileNotFoundError:
        print(f"错误: 找不到映射文件 '{mapping_file}'。请检查文件路径和名称是否正确。")
        sys.exit(1)
    except Exception as e:
        print(f"错误: 读取映射文件时发生错误: {e}")
        sys.exit(1)


def rename_fasta_sequences(fasta_file, output_file, name_map):
    """
    读取FASTA文件，替换序列名称，并生成一个新的单行格式FASTA文件。
    """
    sequences_to_write = []
    replaced_count = 0
    total_sequences = 0

    try:
        # 步骤 1: 解析FASTA文件并替换名称，将结果暂存到列表中
        for record in SeqIO.parse(fasta_file, "fasta"):
            total_sequences += 1
            original_id = record.id
            
            if original_id in name_map:
                new_id = name_map[original_id]
                record.id = new_id
                record.description = "" 
                replaced_count += 1
            
            sequences_to_write.append(record)
        
        # 步骤 2: 手动写入输出文件，以确保每条序列内容只占一行
        with open(output_file, 'w', encoding='utf-8') as f_out:
            for record in sequences_to_write:
                # 写入大于号和新的序列ID，然后换行
                f_out.write(f">{record.id}\n")
                # 写入完整的序列内容，然后换行
                f_out.write(f"{record.seq}\n")

        print(f"成功: 已处理文件 '{fasta_file}'。")
        print(f"统计: 共检查了 {total_sequences} 条序列，成功替换了 {replaced_count} 个名称。")
        print(f"结果已保存到 '{output_file}' (单行格式)。")

    except FileNotFoundError:
        print(f"错误: 找不到FASTA文件 '{fasta_file}'。请检查文件路径和名称是否正确。")
        sys.exit(1)
    except Exception as e:
        print(f"错误: 处理FASTA文件时发生错误: {e}")
        print("请确保你的输入文件是有效的FASTA格式。")
        sys.exit(1)


def main():
    """
    主函数，执行整个流程。
    """
    print("--- 开始批量替换FASTA序列名称并转换为单行格式 ---")
    
    name_replacement_map = create_replacement_map(MAPPING_FILE_PATH)
    
    if name_replacement_map:
        rename_fasta_sequences(FASTA_FILE_PATH, OUTPUT_FASTA_PATH, name_replacement_map)
        
    print("--- 任务完成 ---")


if __name__ == '__main__':
    main()