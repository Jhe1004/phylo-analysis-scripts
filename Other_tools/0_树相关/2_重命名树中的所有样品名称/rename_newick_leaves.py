# -*- coding: utf-8 -*-

import csv
import sys
from ete3 import Tree

# =============================================================================
#                            用户配置区域
#
#  请在这里修改你的输入和输出文件路径。
#  你只需要更改引号内的文件名即可。
# =============================================================================

# 1. 输入的Newick树文件路径 (你的原始树文件)
#    这个文件可以包含一个或多个Newick格式的树，每个树占一行。
TREE_FILE_PATH = 'RAxML_bestTree.result.tree'

# 2. 名称映射CSV文件路径
#    CSV文件格式:
#    第1列: 原始物种名 (在树文件中显示的名字)
#    第2列: 你想要替换成的新物种名
#    注意: CSV文件不应包含表头(header)。
MAPPING_FILE_PATH = 'mapping.csv'

# 3. 输出的Newick树文件路径 (结果将保存到这里)
#    脚本将创建一个新文件，其中包含名称更新后的树。
OUTPUT_FILE_PATH = 'renamed_tree.newick'

# =============================================================================
#                            脚本核心逻辑
#
#            一般来说，你不需要修改下面的代码。
# =============================================================================

def create_replacement_map(mapping_file):
    """
    从CSV文件中读取名称映射关系，并创建一个字典。
    
    参数:
        mapping_file (str): 映射CSV文件的路径。
        
    返回:
        dict: 一个将旧名称映射到新名称的字典。
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
        sys.exit(1) # 退出脚本
    except Exception as e:
        print(f"错误: 读取映射文件时发生错误: {e}")
        sys.exit(1)


def rename_leaves_in_tree(tree_file, output_file, name_map):
    """
    读取Newick文件，替换叶子节点（物种）名称，并写入新文件。
    
    参数:
        tree_file (str): 输入的Newick文件路径。
        output_file (str): 输出的Newick文件路径。
        name_map (dict): 名称映射字典。
    """
    replaced_count = 0
    total_leaves = 0
    
    try:
        with open(tree_file, 'r', encoding='utf-8') as infile, \
             open(output_file, 'w', encoding='utf-8') as outfile:
            
            # 逐行读取，因为一个文件可能包含多个树
            for i, newick_string in enumerate(infile):
                newick_string = newick_string.strip()
                if not newick_string:
                    continue # 跳过空行

                # 使用ete3加载树，格式`1`保留分支长度和内部节点名称
                tree = Tree(newick_string, format=1)
                
                # 遍历树的所有叶子节点
                for leaf in tree.iter_leaves():
                    total_leaves += 1
                    # 如果当前叶子的名称在我们的映射字典中，就替换它
                    if leaf.name in name_map:
                        leaf.name = name_map[leaf.name]
                        replaced_count += 1
                
                # 将修改后的树以Newick格式写回输出文件
                # format=1确保保留所有原始信息（分支长度、节点名等）
                outfile.write(tree.write(format=1))
                outfile.write('\n') # 每个树占一行

        print(f"成功: 已处理文件 '{tree_file}'。")
        print(f"统计: 共检查了 {total_leaves} 个叶子节点，成功替换了 {replaced_count} 个名称。")
        print(f"结果已保存到 '{output_file}'。")

    except FileNotFoundError:
        print(f"错误: 找不到树文件 '{tree_file}'。请检查文件路径和名称是否正确。")
        sys.exit(1)
    except Exception as e:
        print(f"错误: 处理树文件时发生错误: {e}")
        print("请确保你的输入文件是有效的Newick格式。")
        sys.exit(1)

def main():
    """
    主函数，执行整个流程。
    """
    print("--- 开始批量替换Newick树物种名称 ---")
    
    # 步骤1: 创建名称映射字典
    name_replacement_map = create_replacement_map(MAPPING_FILE_PATH)
    
    # 步骤2: 替换并保存新树
    if name_replacement_map:
        rename_leaves_in_tree(TREE_FILE_PATH, OUTPUT_FILE_PATH, name_replacement_map)
        
    print("--- 任务完成 ---")


if __name__ == '__main__':
    main()