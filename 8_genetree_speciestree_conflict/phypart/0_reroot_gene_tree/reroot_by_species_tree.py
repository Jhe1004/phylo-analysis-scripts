import os
from ete3 import Tree
now_dir = os.getcwd()    

species_tree_file = "cds_min_seq_10.tree.newick"  #输入一棵物种树，所有基因树都会根据这棵物种树进行置根
ingroup_file = "ingroups.txt"  #输入想要分析的ingroups物种的名称文件，文件的每一行都表示一个物种


#获取当前文件夹中指定文件
def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if "RAxML_bipartitions" in each:
            file_list.append(each)
    return file_list


def get_outgroup_sort_list():
    species_tree_outgroup_sort_list = []
    species_tree_ingroup_list = []
    with open(ingroup_file, "r") as read_file:
        for each_line in read_file:
            if len(each_line) > 2:
                species_tree_ingroup_list.append(each_line.replace("\n", ""))  
    species_tree = Tree(species_tree_file)
    def get_node(node):
        if node.is_leaf():
            pass
        else:
            if len(species_tree_ingroup_list) == len(node.get_leaf_names()):
                pass
            else:
                if species_tree_ingroup_list[0] in node.children[0].get_leaf_names():
                    species_tree_outgroup_sort_list.append(node.children[1].get_leaf_names())
                    node = node.children[0]
                else:
                    species_tree_outgroup_sort_list.append(node.children[0].get_leaf_names())
                    node = node.children[1]
                get_node(node)
    get_node(species_tree)
    return species_tree_outgroup_sort_list

def get_number_speciestree_outgroups_in_genetree_ingroups(ingroup_node, species_tree_ingroup_list):
    n = 0
    for each_node in ingroup_node.get_leaf_names():
        if each_node not in species_tree_ingroup_list:
            n = n + 1
    return n


def reroot(each_tree):
    print("gene tree: " + each_tree)
    #输入物种树
    species_tree = Tree(species_tree_file)
    #输入基因树
    gene_tree = Tree(each_tree)
    #获取物种树种ingroups的列表
    species_tree_ingroup_list = []
    with open(ingroup_file, "r") as read_file:
        for each_line in read_file:
            if len(each_line) > 2:
                species_tree_ingroup_list.append(each_line.replace("\n", "")) 
    #获取species_tree_outgroup_sort_list  
    species_tree_outgroup_sort_list = get_outgroup_sort_list()
    #获取这棵基因树的ingroup列表和outgroup列表
    gene_tree_group_list = []
    gene_tree_ingroup_list = []
    gene_tree_outgroup_list = []
    for each_leaf in gene_tree.get_leaf_names():
        gene_tree_group_list.append(each_leaf)
        if each_leaf in species_tree_ingroup_list:
            gene_tree_ingroup_list.append(each_leaf)
        else:
            gene_tree_outgroup_list.append(each_leaf)
    #从物种树最外面的根开始遍历
    for each_hypothesis_root in species_tree_outgroup_sort_list:
        if len(each_hypothesis_root) == 1:
            if each_hypothesis_root[0] in gene_tree_group_list:
                gene_tree.set_outgroup(each_hypothesis_root[0])
                ancestor_node = gene_tree.get_common_ancestor(gene_tree_ingroup_list)
                a = get_number_speciestree_outgroups_in_genetree_ingroups(ancestor_node, species_tree_ingroup_list)
                if a/len(gene_tree_outgroup_list) < 0.3333:
                    print("out from aa")
                    return gene_tree
                else:
                    pass
            else:
                pass
        else:
            mono_list = []
            for each in each_hypothesis_root:
                if each in gene_tree_group_list:
                    mono_list.append(each)
            if len(mono_list) == 0:
                pass
            elif len(mono_list) == 1:
                gene_tree.set_outgroup(mono_list[0])
                ancestor_node = gene_tree.get_common_ancestor(gene_tree_ingroup_list)
                a = get_number_speciestree_outgroups_in_genetree_ingroups(ancestor_node, species_tree_ingroup_list)
                if a/len(gene_tree_outgroup_list) < 0.3333:
                    print("out from bb")
                    return gene_tree
                else:
                    pass
            else:
                outgroup_ancestor = gene_tree.get_common_ancestor(mono_list)
                if len(outgroup_ancestor.get_leaf_names()) == len(mono_list):
                    gene_tree.set_outgroup(outgroup_ancestor)
                    ancestor_node = gene_tree.get_common_ancestor(gene_tree_ingroup_list)
                    a = get_number_speciestree_outgroups_in_genetree_ingroups(ancestor_node, species_tree_ingroup_list)
                    if a/len(gene_tree_outgroup_list) < 0.3333:
                        print("out from cc")
                        return gene_tree
                    else:
                        pass
                else:
                    for each in mono_list:
                        gene_tree.set_outgroup(each)
                        ancestor_node = gene_tree.get_common_ancestor(gene_tree_ingroup_list)
                        a = get_number_speciestree_outgroups_in_genetree_ingroups(ancestor_node, species_tree_ingroup_list)
                        if a/len(gene_tree_outgroup_list) < 0.3333:
                            print("out from dd")
                            return gene_tree
                        else:
                            pass

    
      

def main():
    fasta_file_list = get_file_list()
    for each_tree in fasta_file_list:
        if reroot(each_tree):
            reroot(each_tree).write(format=1, outfile=each_tree + "_r.tree")
        
        

if __name__ == "__main__":
    main()