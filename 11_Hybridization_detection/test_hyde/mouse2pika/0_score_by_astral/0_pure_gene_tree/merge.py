import os

def merge_new_tree_files(output_filename="merged_new_trees.txt"):
    # 获取当前目录
    current_dir = os.getcwd()
    
    # 获取所有 *new.tree 文件
    tree_files = [f for f in os.listdir(current_dir) if f.endswith("new.tree") and os.path.isfile(f)]
    
    # 检查是否找到文件
    if not tree_files:
        print("No files matching '*new.tree' found in the current directory.")
        return
    
    # 按文件名排序（可选，确保合并顺序一致）
    tree_files.sort()
    
    # 合并文件
    with open(output_filename, "w") as outfile:
        for i, filename in enumerate(tree_files):
            with open(filename, "r") as infile:
                content = infile.read().rstrip('\n')  # 移除文件末尾的换行符
                outfile.write(content)
            # 在除最后一个文件外的每个文件后添加换行符
            if i < len(tree_files) - 1:
                outfile.write("\n")
    
    print(f"Successfully merged {len(tree_files)} '*new.tree' files into {output_filename}")

if __name__ == "__main__":
    merge_new_tree_files()