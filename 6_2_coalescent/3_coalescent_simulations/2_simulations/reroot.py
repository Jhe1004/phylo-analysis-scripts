'''
使用R语言包phybase中的功能sim.coaltree.sp，基于astral生成的物种树对基因树进行模拟，共模拟10000棵基因树，最后再将这些基因树置根
输入：*.regular.tre格式的超度量的Astral树文件
设置总共模拟的次数，也就是最后模拟出多少棵树
输出：simulate.tree文件，其中包含了所有模拟出的树
'''
import os
from ete3 import Tree


outgroup_species = "60CD.fasta.transdecoder.pep"


#设置外类群
with open("simulate.trees", "r") as read_file:
    with open("astral.simulate.reroot.tree", "a") as write_file:
        for each_line in read_file:
            tree1 = Tree(each_line)
            tree1.set_outgroup(outgroup_species)
            write_file.write(tree1.write() + "\n")


