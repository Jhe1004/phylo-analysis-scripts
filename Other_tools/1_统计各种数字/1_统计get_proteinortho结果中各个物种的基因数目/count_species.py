# count_species.py
import os
from collections import Counter

# --- 配置 ---
# 您要查找的文件后缀名
file_suffix = 'cds.fasta'
# --- 配置结束 ---

def count_species_in_fasta_files(path='.', suffix='.fasta'):
    """
    统计指定路径下具有特定后缀的FASTA文件中各物种的出现频率。

    Args:
        path (str): 要搜索的文件夹路径，默认为当前文件夹。
        suffix (str): FASTA文件的后缀名。

    Returns:
        collections.Counter: 一个计数器对象，其中键是物种名，值是该物种出现的（不同）文件数。
    """
    # 1. 查找所有符合条件的fasta文件
    try:
        fasta_files = [f for f in os.listdir(path) if f.endswith(suffix)]
    except FileNotFoundError:
        print(f"错误：找不到文件夹 '{path}'。请确保脚本在正确的路径下运行。")
        return None

    if not fasta_files:
        print(f"在当前文件夹中没有找到任何以 '{suffix}' 结尾的文件。")
        return None

    print(f"成功找到 {len(fasta_files)} 个文件进行分析...")

    # 2. 初始化一个总的计数器
    # 使用Counter可以非常方便地进行计数
    total_species_counts = Counter()

    # 3. 遍历每一个找到的fasta文件
    for filename in fasta_files:
        # 对于每个文件，我们只希望对一个物种计数一次，即使它在该文件中有多个序列。
        # 使用集合(set)可以自动处理重复项。
        species_in_current_file = set()
        
        try:
            with open(os.path.join(path, filename), 'r', encoding='utf-8') as f:
                for line in f:
                    # FASTA格式的序列标识符以 '>' 开头
                    if line.startswith('>'):
                        # 假设物种名是'>'之后到第一个空格或换行符之间的字符串
                        # 例如：从 ">Homo_sapiens gene_xyz" 中提取 "Homo_sapiens"
                        header = line.strip()
                        # [1:] 去掉'>', .split()[0] 取第一部分
                        species_name = header.split()[0][1:]
                        species_in_current_file.add(species_name)
        except Exception as e:
            print(f"读取文件 {filename} 时发生错误: {e}")
            continue
            
        # 4. 将当前文件中出现过的所有（不重复的）物种，更新到总计数器中
        # 这一步等价于 for species in species_in_current_file: total_species_counts[species] += 1
        total_species_counts.update(species_in_current_file)

    return total_species_counts

def print_results(species_counts):
    """
    格式化并打印结果。
    """
    if not species_counts:
        print("\n分析完成，但没有找到任何物种信息。")
        print("请检查您的FASTA文件是否为空，或者序列标识符（以'>'开头的行）格式是否正确。")
        return

    print("\n--- 统计结果 ---")
    print(f"{'物种 (Species)':<30} {'出现过的文件数 (Count)':<20}")
    print("-" * 50)

    # 按照出现次数从高到低排序后输出
    for species, count in species_counts.most_common():
        print(f"{species:<30} {count:<20}")

# --- 脚本主程序 ---
if __name__ == "__main__":
    # 执行主函数
    final_counts = count_species_in_fasta_files(path='.', suffix=file_suffix)
    
    # 如果分析成功，则打印结果
    if final_counts is not None:
        print_results(final_counts)