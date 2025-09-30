#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from ete3 import Tree

# --- 参数配置区域 ---
# 您可以在这里修改所有参数，无需在命令行输入

# 1. 定义要搜索的目标文件的后缀
TARGET_SUFFIX = ".trees"

# 2. 定义输出文件夹的名称
OUTPUT_DIR = "rerooted_trees"

# 3. 定义外群（Outgroup）物种列表
OUTGROUP_LIST = [
    "Clematis_gratopsis.fasta.transdecoder.pep"
]

# --- 函数定义 ---

def get_target_files():
    """
    获取当前目录下所有以 TARGET_SUFFIX 结尾的文件名。
    """
    print(f"正在搜索以 '{TARGET_SUFFIX}' 结尾的文件...")
    try:
        # 使用在配置区定义的 TARGET_SUFFIX 变量
        file_list = [f for f in os.listdir('.') if f.endswith(TARGET_SUFFIX)]
        return file_list
    except OSError as e:
        print(f"错误：无法读取当前目录下的文件。 {e}", file=sys.stderr)
        return []

def reroot_tree_file(input_filename, output_dir):
    """
    读取指定的树文件，根据外群列表进行置根，并将成功置根的树写入到输出目录的新文件中。

    Args:
        input_filename (str): 原始树文件的名称。
        output_dir (str): 输出文件夹的路径。
    """
    print(f"--- 开始处理文件: {input_filename} ---")
    
    # 构建输出文件的完整路径
    output_filepath = os.path.join(output_dir, os.path.basename(input_filename) + "_reroot.tree")
    
    successful_reroot_count = 0
    
    try:
        with open(input_filename, "r") as read_file:
            for line_num, each_line in enumerate(read_file, 1):
                line_content = each_line.strip()
                if not line_content:
                    continue

                try:
                    t = Tree(line_content)
                    outgroup_nodes = t.get_leaves_by_name(OUTGROUP_LIST)
                    
                    if outgroup_nodes:
                        ancestor = t.get_common_ancestor(outgroup_nodes)
                        t.set_outgroup(ancestor)
                        
                        with open(output_filepath, "a") as write_file:
                            write_file.write(t.write() + "\n")
                        
                        successful_reroot_count += 1
                    else:
                        pass

                except Exception as e:
                    print(f"  > 警告: 处理文件 '{input_filename}' 的第 {line_num} 行时出错: {e}")

    except IOError as e:
        print(f"错误：无法读取文件 '{input_filename}'。 {e}", file=sys.stderr)
        return

    if successful_reroot_count > 0:
        print(f"✅ 处理完成: {successful_reroot_count} 棵树被成功置根并保存到 '{output_filepath}'")
    else:
        print(f"❌ 处理完成: 在 '{input_filename}' 中没有找到任何可以被成功置根的树。")


def main():
    """
    脚本的主执行流程。
    """
    # 检查并安装依赖
    try:
        from ete3 import Tree
    except ImportError:
        print("错误：核心依赖 'ete3' 模块未安装。", file=sys.stderr)
        print("请使用命令安装: pip install ete3", file=sys.stderr)
        sys.exit(1)
        
    # 1. 创建输出文件夹
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"已创建输出文件夹: '{OUTPUT_DIR}'")

    # 2. 获取要处理的文件列表
    tree_file_list = get_target_files()
    if not tree_file_list:
        print(f"\n在当前目录中没有找到任何以 '{TARGET_SUFFIX}' 结尾的文件。")
        return

    print(f"\n找到 {len(tree_file_list)} 个目标文件，开始进行置根处理...")
    
    # 3. 循环处理每个文件
    for each_file in tree_file_list:
        reroot_tree_file(each_file, OUTPUT_DIR)
    
    print("\n---------------------------------")
    print("所有文件处理完毕！")
    print("---------------------------------")


# 脚本入口
if __name__ == "__main__":
    main()