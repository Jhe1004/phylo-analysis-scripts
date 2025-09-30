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

input_tree_file = "result.tree"
output_tree_file = "out.trees"
subtree_species_file = "list.txt"
outgroup_species = "Gossypium_hirsutum_SRR8156069.fasta.transdecoder.pep"


def main():
    subtree_species_list = []
    with open(subtree_species_file, "r") as read_file:
        for each_line in read_file:
            if len(each_line) > 2:
                subtree_species_list.append(each_line.replace("\n", "").replace(" ", ""))

    with open(output_tree_file, "a") as write_file:
        with open(input_tree_file, "r") as read_file:
            for each_tree in read_file:
                tree1 = Tree(each_tree)
                try:
                    tree1.prune(subtree_species_list)
                    tree1.set_outgroup(outgroup_species)
                    write_file.write(tree1.write() + "\n")
                except :
                    continue

main()