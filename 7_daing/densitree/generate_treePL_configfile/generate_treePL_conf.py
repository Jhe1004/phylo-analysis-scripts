'''
批量的生成treePL软件所需要的配置文件
输入：.tree后缀的树文件，个文件中含有一棵树，所有树需要含有完全相同的物种名称
需要指定根节点的一个精确时间，一般就是使用物种树dating完之后根节点的精确时间
输出：treePL的配置文件
'''
import os
from ete3 import Tree
import dendropy

now_dir = os.getcwd()    


calibrated_tree = "dating.tre"

#函数get_file_list：获取当前文件夹中指定文件
def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if ".new" in each:
            file_list.append(each)
    return file_list


def get_calibrat_age(monophy):
    tree1 =  dendropy.Tree.get_from_path(calibrated_tree, "newick", preserve_underscores=True)
    tree1.calc_node_ages()
    mrca = tree1.mrca(taxon_labels=monophy)
    return(mrca.age)



def generate_calibra_list(each_tree):
    calibra_list = []
    t = Tree(each_tree)
    calibrated = Tree(calibrated_tree)
    def check_monophyly(t, calibrated, each_node):
        gene_tree_ingroup = each_node.get_leaf_names()
        gene_tree_outgroup = []
        for each in t.get_leaf_names():
            if each not in gene_tree_ingroup:
                gene_tree_outgroup.append(each)
        for each_leaf in calibrated.get_common_ancestor(gene_tree_ingroup).get_leaf_names():
            if each_leaf in gene_tree_outgroup:
                return False
        return True      
    #calibra root
    calibra_list.append("mrca = a1 " + str(t.get_leaf_names()).replace("\'", "").replace(",", "")[1:-1] + "\n")
    calibra_list.append("min = a1 " + str(get_calibrat_age(t.get_leaf_names())) + "\n")
    calibra_list.append("max = a1 " + str(get_calibrat_age(t.get_leaf_names()) + 0.000001) + "\n")
    #calibra others
    n = 27
    for each_node in t.traverse():
        if each_node.is_leaf():
            pass
        else:
            if n == 28:
                n = n + 1
            elif check_monophyly(t, calibrated, each_node):
                calibra_list.append("mrca = a" + str(n) + " " + str(each_node.get_leaf_names()).replace("\'", "").replace(",", "")[1:-1] + "\n")
                calibra_list.append("min = a" + str(n) + " " + str(get_calibrat_age(each_node.get_leaf_names()) - 0.0005) + "\n")
                calibra_list.append("max = a" + str(n) + " " + str(get_calibrat_age(each_node.get_leaf_names()) + 0.005) + "\n")
                n = n + 1
    return calibra_list





#生成treePL配置文件
def generate_treePL_conf(each_tree):
    with open(each_tree.replace("RAxML_bipartitions.", "")[:-4] + ".fasta") as read_file:
        for each_line in read_file:
            if each_line[0] != ">":
                seq_len = len(each_line)-1
                break
    with open(each_tree[:-5] + ".conf", "a") as write_file:
        write_file.write("[Input files containing the ML trees]\n")
        write_file.write("treefile = " + each_tree + "\n")         
        write_file.write("[General commands]\n"
                        "nthreads = 12\n"
                        "smooth = 0.0000001\n"
                        "thorough\n"
                        "log_pen\n")
        write_file.write("numsites = "  + str(seq_len) + "\n")

        write_file.write("[Calibrations]\n")
        for each_line in generate_calibra_list(each_tree):
            write_file.write(each_line)

        write_file.write("[Best optimisation parameters]\n")
        write_file.write("opt = 3\n"
                        "moredetail\n"
                        "optad = 3\n"
                        "moredetailad\n"
                        "optcvad = 5\n"          
                        "[randomcv]\n"
                        "[cviter = 5]\n"
                        "[cvsimaniter = 100000]\n"
                        "[cvstart = 10000]\n"
                        "[cvstop = 0.00001]\n"
                        "[cvmultstep = 0.1]\n"
                        "cvoutfile = " + each_tree + ".cvout\n"
                        "outfile = " + each_tree + ".newtree\n")         

    

def main():
    fasta_file_list = get_file_list()
    for each_tree in fasta_file_list:
        generate_treePL_conf(each_tree)
        

if __name__ == "__main__":
    main()