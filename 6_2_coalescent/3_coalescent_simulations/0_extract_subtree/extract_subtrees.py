import os
from ete3 import Tree
import argparse


parser = argparse.ArgumentParser(description="", add_help=True)
required = parser.add_argument_group("Required arguments")
required.add_argument('-i', '--input_gene_file', action="store", metavar='\b', type=str, required=True, help="")  
required.add_argument('-n', '--subtree_species_file', action="store", metavar='\b', type=str, required=True, help="")
required.add_argument('-o', '--outgroup_species', action="store", metavar='\b', type=str, required=True, help="")

args                           = parser.parse_args()
input_gene_file                = args.input_gene_file
subtree_species_file           = args.subtree_species_file
outgroup_species               = args.outgroup_species


subtree_species_list = []
with open(subtree_species_file, "r") as read_file:
    for each_line in read_file:
        if len(each_line) > 2:
            subtree_species_list.append(each_line.replace("\n", ""))

tree_list = []
with open(input_gene_file, "r") as read_file:
    for each_line in read_file:
        if len(each_line) > 2:
            tree_list.append(each_line.replace("\n", ""))

with open("output_subtree.trees", "a") as write_file:
    for each_tree in tree_list:
        tree1 = Tree(each_tree)
        try:
            tree1.prune(subtree_species_list)
            tree1.set_outgroup(outgroup_species)
            write_file.write(tree1.write() + "\n")
        except:
            continue

