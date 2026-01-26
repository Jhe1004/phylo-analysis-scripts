'''
比较几棵基因树两两之间的距离
输入的是已经置根的树，后缀为newick
输出一个距离表格
'''
import os
import numpy as np
import pandas as pd
import dendropy
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt
import seaborn as sns
now_dir = os.getcwd() 


gene_trees = "output_subtree.trees"
simulated_trees = "astral.simulate.reroot.tree"
astral_tree = "astral.tre.regular.tre"


def get_distance(tree1, tree2):
    tns = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get_from_string(tree1, "newick", taxon_namespace=tns)
    tree2 = dendropy.Tree.get_from_string(tree2, "newick", taxon_namespace=tns)
    return treecompare.symmetric_difference(tree1, tree2)

def compare(trees1, astral_tree):
    with open(astral_tree, "r") as read_file:
        astral_tree = read_file.readline().replace("\n", "")
    compare_list = []
    with open(trees1, "r") as read_file:
        for each_line in read_file:
            if len(each_line) > 2:
                compare_list.append(get_distance(each_line.replace("\n", ""), astral_tree))
    compare_prob_density_dict = {}
    compare_list.sort()
    for each_num in range(0, compare_list[-1] + 2, 2):
        compare_prob_density_dict[each_num] = compare_list.count(each_num)/len(compare_list)
    return compare_prob_density_dict
            

def main():
    genetree_prob_density = compare(gene_trees, astral_tree)
    simulatedtree_prob_density = compare(simulated_trees, astral_tree)
    with open("result.tsv", "w") as write_file:
        write_file.write("distance\tclass\tprob\n")
        for each in genetree_prob_density:
            write_file.write(str(each) + "\t" + gene_trees + "\t" + str(genetree_prob_density[each]) + "\n")
        for each in simulatedtree_prob_density:
            write_file.write(str(each) + "\t" + simulated_trees + "\t" + str(simulatedtree_prob_density[each]) + "\n")
    compare_data = pd.read_csv("result.tsv", sep="\t")
    f, ax = plt.subplots(figsize=(20, 10))
    sns.barplot(x="distance", y="prob", hue="class", data=compare_data)
    plt.savefig("compare.png", dpi = 120)
    




main()