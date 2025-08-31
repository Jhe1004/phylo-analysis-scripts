import sys
from Bio import AlignIO
from collections import Counter

def calculate_alignment_stats(file_path):
    """
    计算FASTA比对文件的统计信息，包括总长度、变异位点和简约信息位点。

    参数:
    file_path (str): FASTA格式比对文件的路径。

    返回:
    dict: 包含统计信息的字典。
    """
    try:
        # 读取比对文件
        alignment = AlignIO.read(file_path, "fasta")
    except FileNotFoundError:
        print(f"错误：找不到文件 '{file_path}'")
        return None
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return None

    # 获取比对长度（所有序列的长度应该是一样的）
    alignment_length = alignment.get_alignment_length()

    # 初始化计数器
    variable_sites = 0
    parsimony_informative_sites = 0

    # 逐个位点（列）进行分析
    for i in range(alignment_length):
        # 获取当前位置的所有字符
        column = alignment[:, i]
        
        # 移除gap '-' 和未知字符 'N' or '?' 等，只统计ATCG
        # 如果您想将gap也视为一种变异，可以注释掉下面这行
        column = ''.join([base for base in column.upper() if base in 'ATCG'])

        if not column:
            continue # 如果该列只有gap或未知字符，则跳过

        # 计算该位点的字符频率
        counts = Counter(column)
        unique_bases = [base for base in counts if base in 'ATCG']

        # 检查是否是变异位点
        # 如果存在至少两种不同的碱基，则为变异位点
        if len(unique_bases) > 1:
            variable_sites += 1

            # 检查是否是简约信息位点 (Parsimony-Informative Site)
            # 条件：至少有两种字符，且每种字符至少出现两次
            if len(unique_bases) >= 2:
                count_of_bases_appearing_at_least_twice = 0
                for base in unique_bases:
                    if counts[base] >= 2:
                        count_of_bases_appearing_at_least_twice += 1
                
                if count_of_bases_appearing_at_least_twice >= 2:
                    parsimony_informative_sites += 1

    return {
        "characters": alignment_length,
        "variable_sites": variable_sites,
        "parsimony_informative_sites": parsimony_informative_sites
    }

if __name__ == "__main__":
    # 检查命令行是否提供了文件名
    if len(sys.argv) != 2:
        print("使用方法: python alignment_stats.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    stats = calculate_alignment_stats(fasta_file)

    if stats:
        print("比对文件统计结果:")
        print(f"  总长度 (Characters): {stats['characters']}")
        print(f"  变异位点 (Variable sites): {stats['variable_sites']}")
        print(f"  简约信息位点 (Parsimony-informative sites): {stats['parsimony_informative_sites']}")