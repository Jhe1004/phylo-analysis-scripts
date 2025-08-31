import os
import sys
from ete3 import Tree

# --- 参数配置 ---

# 定义输出文件夹的名称
OUTPUT_DIR = "rerooted_trees"

# 定义外群（Outgroup）物种列表
# 您可以根据需要在此列表中添加或修改物种名称
OUTGROUP_LIST = [
    "Clematis_gratopsis.fasta.transdecoder.pep"
]

# --- 函数定义 ---

def get_target_files():
    """
    获取当前目录下所有以 '_cds_maffted' 结尾的文件名。
    """
    try:
        # 使用列表推导式更简洁
        file_list = [f for f in os.listdir('.') if f.endswith("_cds_maffted")]
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
    output_filepath = os.path.join(output_dir, input_filename + "_reroot.tree")
    
    successful_reroot_count = 0
    
    try:
        with open(input_filename, "r") as read_file:
            # 逐行读取（即逐个处理树）
            for line_num, each_line in enumerate(read_file, 1):
                line_content = each_line.strip()
                
                # 跳过空行
                if not line_content:
                    continue

                try:
                    t = Tree(line_content)

                    # 查找外群节点
                    # get_leaves_by_name可以接受单个名称或名称列表
                    outgroup_nodes = t.get_leaves_by_name(OUTGROUP_LIST)
                    
                    if outgroup_nodes:
                        # 如果找到多个外群节点，使用它们的最近公共祖先(LCA)
                        # 如果只找到一个，LCA就是它本身
                        ancestor = t.get_common_ancestor(outgroup_nodes)
                        t.set_outgroup(ancestor)
                        
                        # 只有在成功置根后，才以追加模式打开并写入文件
                        with open(output_filepath, "a") as write_file:
                            write_file.write(t.write() + "\n")
                        
                        successful_reroot_count += 1
                    else:
                        # 在文件中未找到任何指定的外群
                        # print(f"  > 在第 {line_num} 行的树中未找到外群，跳过。")
                        pass

                except Exception as e:
                    print(f"  > 警告: 处理文件 '{input_filename}' 的第 {line_num} 行时出错: {e}")

    except IOError as e:
        print(f"错误：无法读取文件 '{input_filename}'。 {e}", file=sys.stderr)
        return

    # 根据处理结果打印最终信息
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
        print("\n在当前目录中没有找到任何以 '_cds_maffted' 结尾的文件。")
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
