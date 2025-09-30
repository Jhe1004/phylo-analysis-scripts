import os
from ete3 import Tree
import dendropy

now_dir = os.getcwd()


calibrated_tree = "calibrated_tree"
calibrated_species = "Berberis_thunbergii.fasta.transdecoder.pep"
ingroupspecies_list_file = "ingroupspecies.txt"


#函数get_file_list：获取当前文件夹中指定文件
def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:
        if ".tree" in each:
            file_list.append(each)
    return file_list

#获取基因树一个节点后的所有物种名称
def get_species_in_clade(each_node):
    species_in_clade_list = []
    for each in each_node.get_leaf_names():
        if each.split("++")[0] not in species_in_clade_list:
            species_in_clade_list.append(each.split("++")[0])
    return species_in_clade_list


#根据一些树中的物种名称，获得标定时间
def get_calibrat_age(monophy):
    t = Tree(calibrated_tree)
    b = t.get_common_ancestor(monophy)
    c = b.get_leaf_names()

    tree1 =  dendropy.Tree.get_from_path(calibrated_tree, "newick", preserve_underscores=True)
    tree1.calc_node_ages()
    mrca = tree1.mrca(taxon_labels=c)
    return(mrca.age)


#生成标定列表
def generate_calibra_list(each_tree):
    calibra_list = []
    t = Tree(each_tree)
    #获取ingroups_list
    ingroups_list = []
    with open(ingroupspecies_list_file, "r") as read_file:
        for each_ingroup_species in read_file:
            ingroups_list.append(each_ingroup_species.replace("\n",""))

    #首先判断Berberis在不在树中，不在则return None
    if calibrated_species in str(t.get_leaf_names()): #berberis 在树中
        for each_node in t.traverse():
            if each_node.children: #如果不是叶子节点
                if get_species_in_clade(each_node.children[0])[0] == calibrated_species and len(get_species_in_clade(each_node.children[0])) == 1: #如果该节点是berberis的父节点
                    # 判断berberis的姊妹群中的物种是否都是ingroups里面的
                    judge = True
                    for each in get_species_in_clade(each_node.children[1]):
                        if each not in ingroups_list:
                            judge = False
                    if judge:
                        # 判断该姊妹群中包含的物种数是否占整个基因树内类群物种的1/2以上(为了保证berberis的系统位置基本正确)
                        a = []
                        for each in  get_species_in_clade(each_node.children[1]):
                            if each in ingroups_list:
                                a.append(each)
                        b = []
                        for each in  get_species_in_clade(t):
                            if each in ingroups_list:
                                b.append(each)
                        if len(a)/len(b) > 2/3:
                            # 再判断该节点是不是root，如果是root则标定根节点，不是root则标定该节点以及根节点
                            calibra_list.append("mrca = a1 " + str(t.get_leaf_names()).replace("\'", "").replace(",", "")[1:-1] + "\n")
                            calibra_list.append("min = a1 " + str(get_calibrat_age(get_species_in_clade(t))) + "\n")
                            calibra_list.append("max = a1 " + str(get_calibrat_age(get_species_in_clade(t)) + 0.000001) + "\n")   
                            if each_node.is_root():
                                pass
                            else:
                                calibra_list.append("mrca = a2 " + str(each_node.get_leaf_names()).replace("\'", "").replace(",", "")[1:-1] + "\n")
                                calibra_list.append("min = a2 " + str(get_calibrat_age(get_species_in_clade(each_node))) + "\n")
                                calibra_list.append("max = a2 " + str(get_calibrat_age(get_species_in_clade(each_node)) + 0.000001) + "\n")                                
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                pass
    else:
        pass
    return calibra_list


#生成treePL配置文件
def generate_treePL_conf(each_tree):
    calibra_list = generate_calibra_list(each_tree)
    if calibra_list:
        with open(each_tree.replace("RAxML_bipartitions.", "")[:-7] + ".fasta") as read_file:
            for each_line in read_file:
                if each_line[0] != ">":
                    seq_len = len(each_line)-1
                    break
        with open(each_tree[:-5] + ".conf", "a") as write_file:
            write_file.write("[Input files containing the ML trees]\n")
            write_file.write("treefile = " + each_tree + "\n")         
            write_file.write("[General commands]\n"
                            "nthreads = 12\n"
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
                            "smooth = 0.0000001\n"        
                            "[randomcv]\n"
                            "[cviter = 5]\n"
                            "[cvsimaniter = 10000]\n"
                            "[cvstart = 1000]\n"
                            "[cvstop = 0.000001]\n"
                            "[cvmultstep = 0.1]\n"
                            "cvoutfile = " + each_tree + ".cvout\n"
                            "outfile = " + each_tree + ".newtree\n")         
                    
                    

            


def main():
    fasta_file_list = get_file_list()
    for each_tree in fasta_file_list:
        generate_treePL_conf(each_tree)


main()