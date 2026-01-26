# -*- coding: utf-8 -*-

import os
from ete3 import Tree
import random

# ==============================================================================
# 参数配置 (在这里修改你的设置)
# ==============================================================================

# 1. 包含物种名称的输入文件名 (每行一个物种)
INPUT_TAXA_FILE = 't.txt'

# 2. 指定用作外类群的物种名称 (必须与输入文件中的某个名称完全匹配)
OUTGROUP_TAXON_NAME = 'Clematis_songorica_Altay'

# 3. 要生成的随机树的数量
NUM_TREES_TO_GENERATE = 10

# 4. 保存输出 Newick 文件的目录名称
OUTPUT_DIRECTORY = 'random_trees'

# ==============================================================================
# 脚本主逻辑 (通常不需要修改以下部分)
# ==============================================================================

def generate_random_trees_with_ete3():
    """
    主函数，使用 ete3 读取物种、生成并保存随机系统发育树。
    """
    print("脚本开始执行 (使用 ete3 库)...")

    # --- 1. 检查并创建输出目录 ---
    try:
        os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)
        print(f"输出目录 '{OUTPUT_DIRECTORY}' 已准备好。")
    except OSError as e:
        print(f"错误：无法创建目录 '{OUTPUT_DIRECTORY}': {e}")
        return

    # --- 2. 读取物种列表 ---
    try:
        with open(INPUT_TAXA_FILE, 'r', encoding='utf-8') as f:
            taxon_names = [line.strip() for line in f if line.strip()]
        
        if not taxon_names:
            print(f"错误：输入文件 '{INPUT_TAXA_FILE}' 为空或未找到。")
            return
            
        print(f"成功从 '{INPUT_TAXA_FILE}' 文件中读取了 {len(taxon_names)} 个物种。")

    except FileNotFoundError:
        print(f"错误：找不到输入文件 '{INPUT_TAXA_FILE}'。请确保该文件与脚本在同一目录下。")
        return

    # --- 3. 验证外类群是否存在 ---
    if OUTGROUP_TAXON_NAME not in taxon_names:
        print(f"错误：指定的外类群 '{OUTGROUP_TAXON_NAME}' 不在物种列表中。")
        print("请检查 'OUTGROUP_TAXON_NAME' 的拼写是否与列表文件中的完全一致。")
        return

    print("\n开始生成系统发育树...")
    # --- 4. 循环生成并保存树 ---
    for i in range(NUM_TREES_TO_GENERATE):
        # 创建一棵空的树
        t = Tree()
        
        # 使用随机拓扑结构填充树的所有叶子节点
        # random.sample 用于确保每次循环的物种添加顺序是随机的，从而产生不同拓扑
        t.populate(len(taxon_names), names_library=random.sample(taxon_names, len(taxon_names)), random_branches=True)

        # 确保树是二叉树（解决可能出现的多分叉）
        t.resolve_polytomy()

        # 使用指定的外类群来确定树的根
        # get_leaves_by_name 返回一个列表，我们取第一个元素
        try:
            outgroup_leaf = t.get_leaves_by_name(OUTGROUP_TAXON_NAME)[0]
            t.set_outgroup(outgroup_leaf)
        except IndexError:
            # 这个错误理论上不会发生，因为前面已经检查过外类群存在
            print(f"严重错误：在树中找不到外类群 '{OUTGROUP_TAXON_NAME}'。")
            continue

        # 构建输出文件路径
        output_filename = f'random_tree_{i+1}.nwk'
        output_path = os.path.join(OUTPUT_DIRECTORY, output_filename)

        # 以 Newick 格式将树写入文件
        # format=9 表示只输出拓扑结构，不包含枝长和内部节点名称
        t.write(outfile=output_path, format=9)
        
        print(f"已生成并保存树到: {output_path}")

    print(f"\n任务完成！所有 {NUM_TREES_TO_GENERATE} 棵树都已保存在 '{OUTPUT_DIRECTORY}' 文件夹中。")


if __name__ == '__main__':
    generate_random_trees_with_ete3()