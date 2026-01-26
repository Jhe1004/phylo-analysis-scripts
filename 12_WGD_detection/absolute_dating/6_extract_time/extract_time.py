import os
import re
from ete3 import Tree
import dendropy
now_dir = os.getcwd()


#函数get_file_list：获取当前文件夹中指定文件
def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if ".newtree" in each:
            file_list.append(each)
    return file_list


def get_species_in_clade(each_node):
    species_in_clade_list = []
    for each in each_node.get_leaf_names():
        if each.split("++")[0] not in species_in_clade_list:
            species_in_clade_list.append(each.split("++")[0])
    return species_in_clade_list


def get_calibrat_age(monophy, each_tree):
    tree1 =  dendropy.Tree.get_from_path(each_tree, "newick", preserve_underscores=True)
    tree1.calc_node_ages() 
    mrca = tree1.mrca(taxon_labels=monophy)
    return(mrca.age)


def extract_time(each_tree):
    t = Tree(each_tree)
    species_list = get_species_in_clade(t)
    for each_species in species_list:
        #print(each_species)
        for each_node in t.traverse():
            if each_node.is_leaf():
                pass
            else:
                if each_species in str(each_node.children[0].get_leaf_names()) and each_species in str(each_node.children[1].get_leaf_names()):
                    with open(each_species, "a") as write_file:
                        write_file.write(str(get_calibrat_age(each_node.get_leaf_names(), each_tree)) + "\n")





def main():
    fasta_file_list = get_file_list()
    for each_tree in fasta_file_list:
        extract_time(each_tree)


if __name__ == "__main__":
    main()