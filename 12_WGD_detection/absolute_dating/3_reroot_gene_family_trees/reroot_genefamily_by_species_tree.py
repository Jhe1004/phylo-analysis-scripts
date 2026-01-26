import os
from ete3 import Tree
now_dir = os.getcwd()    

species_tree_file = "RAxML_bipartitions.result.newick"  #输入一棵物种树，所有基因树都会根据这棵物种树进行置根
ingroup_file = "ingroupspecies.txt"  #输入想要分析的ingroups物种的名称文件，文件的每一行都表示一个物种


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


#输入一个包含多个物种名的列表，在基因树中获取这些物种的所有基因
def get_genes_in_species(each_hypothesis_root, gene_tree):
    hypothesis_root_genes = []
    for each_gene_name in gene_tree.get_leaf_names():
        for each_species_name in each_hypothesis_root:
            if each_species_name == each_gene_name.split("++")[0]:
                hypothesis_root_genes.append(each_gene_name)
    return hypothesis_root_genes



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

    gene_tree_group_list = []
    gene_tree_ingroup_list = []
    gene_tree_outgroup_list = []
    for each_leaf in gene_tree.get_leaf_names():
        gene_tree_group_list.append(each_leaf)
        if each_leaf.split("++")[0] in species_tree_ingroup_list:
            gene_tree_ingroup_list.append(each_leaf)
        else:
            gene_tree_outgroup_list.append(each_leaf)

    if gene_tree_ingroup_list:
        #获取species_tree_outgroup_sort_list  
        species_tree_outgroup_sort_list = get_outgroup_sort_list()
        #开始置根
        for each_hypothesis_root in species_tree_outgroup_sort_list:
            each_hypothesis_root_genes_list = get_genes_in_species(each_hypothesis_root, gene_tree)
            if each_hypothesis_root_genes_list:
                if len(each_hypothesis_root_genes_list) == 1:
                    gene_tree.set_outgroup(each_hypothesis_root_genes_list[0])
                    return(gene_tree)
                else:
                    gene_tree.set_outgroup(gene_tree_ingroup_list[0]) #先随便使用内类群置个根
                    #检查each_hypothesis_root_genes_list是不是个单系
                    ancestor_node = gene_tree.get_common_ancestor(each_hypothesis_root_genes_list)
                    if len(ancestor_node.get_leaf_names()) == len(each_hypothesis_root_genes_list): #是单系
                        gene_tree.set_outgroup(ancestor_node)
                        return(gene_tree)
                    else:   #不是单系
                        #分别使用each_hypothesis_root_genes_list里面的物种尝试置根，看哪个置根之后比较合适
                        root_gene = each_hypothesis_root_genes_list[0]
                        gene_tree.set_outgroup(each_hypothesis_root_genes_list[0])
                        ancestor_node2 = gene_tree.get_common_ancestor(gene_tree_ingroup_list)
                        ingroups_len = len(ancestor_node2.get_leaf_names())
                        for each_hypothesis_root_gene in each_hypothesis_root_genes_list:
                            gene_tree.set_outgroup(each_hypothesis_root_gene)
                            ancestor_node3 = gene_tree.get_common_ancestor(gene_tree_ingroup_list)
                            ingroups_len2 = len(ancestor_node3.get_leaf_names())
                            if ingroups_len2 < ingroups_len:
                                ingroups_len = ingroups_len2
                                root_gene = each_hypothesis_root_gene
                        gene_tree.set_outgroup(root_gene)
                        return(gene_tree)


def main():
    fasta_file_list = get_file_list()
    for each_tree in fasta_file_list:
        if reroot(each_tree):
            reroot(each_tree).write(format=1, outfile=each_tree + "_r.tree")
        
       

if __name__ == "__main__":
    main()