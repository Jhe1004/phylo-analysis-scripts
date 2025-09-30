from ete3 import Tree
import os
import glob

# 获取当前目录
now_dir = os.getcwd()

def main():
    # 获取当前目录下所有后缀为.tree的文件
    tree_files = glob.glob(os.path.join(now_dir, "*.tree"))

    # 遍历每个.tree文件
    for input_tree in tree_files:
        # 读取树文件
        t = Tree(input_tree)
        
        # 为每个树文件生成对应的.txt文件
        output_txt = input_tree + ".txt"
        with open(output_txt, "w") as write_file:
            for each_line in t.get_leaf_names():
                write_file.write(each_line + "\n")

if __name__ == "__main__":
    main()