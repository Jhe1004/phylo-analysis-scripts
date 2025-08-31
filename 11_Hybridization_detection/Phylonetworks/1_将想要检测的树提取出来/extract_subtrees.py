'''
每棵基因树都提取出一棵特定的subtree， 这棵树必须包含所有预先给出的物种。提取之后对这些树进行重新置根
输入：RaxML输出的树文件，每个树文件中包含一棵基因树，所有树文件均放置在与本脚本相同的文件夹中。
.txt文件，其中每一行代表一个需要提取的物种名称，必须和基因树中的严格对应以及需要
指定一个外类群名称，这里为了方便只能选择一个物种作为外类群。
输出：subtree.trees，包含所有处理好的文件
'''
import os
from ete3 import Tree
now_dir = os.getcwd()    

subtree_species_file = "spcies.txt"
outgroup_species = "E2021"


def get_file_list():
    #函数get_file_list：获取当前文件夹中指定文件
    #输入：无
    #输出：指定文件列表：file_list
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if "RAxML" in each:
            file_list.append(each)
    return file_list


def main():
    subtree_species_list = []
    with open(subtree_species_file, "r") as read_file:
        for each_line in read_file:
            if len(each_line) > 2:
                subtree_species_list.append(each_line.replace("\n", ""))
    tree_list = get_file_list()
    
    for each_tree_file in tree_list:
        
        with open(each_tree_file + ".new", "a") as write_file:
            with open(each_tree_file, "r") as read_file:
                for each_tree in read_file:
                    tree1 = Tree(each_tree)
                    try:
                        tree1.prune(subtree_species_list)
                        #tree1.set_outgroup(outgroup_species)
                        write_file.write(tree1.write() + "\n")
                    except :
                        continue


main()