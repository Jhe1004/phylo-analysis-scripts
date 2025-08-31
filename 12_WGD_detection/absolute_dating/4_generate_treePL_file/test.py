'''
生成treepl软件所需的“.conf”文件。
'''
import os
from ete3 import Tree
import dendropy

calibrated_age = 92
ingroups_file = "ingroupspecies.txt"
outgroups_file = "outgroupspecies.txt"

def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:
        if ".tree" in each:
            file_list.append(each)
    return file_list

#生成标定列表
def generate_calibra_list(each_tree):
    #获取用户设定的内类群列表
    ingroups_list = []
    with open(ingroups_file, "r") as read_file:
        for each_line in read_file:
            ingroups_list.append(each_line.replace("\n",""))
    #获取用户设定的外类群列表
    outgroups_list = []
    with open(outgroups_file, "r") as read_file:
        for each_line in read_file:
            outgroups_list.append(each_line.replace("\n",""))  
    '''
    开始遍历物种树, 满足如下条件的节点将会被标定:
    该节点后面的两个子分支,
    其中一个分支中的类群都属于outgroups_list, 
    另一个分支中的类群都属于ingroups_list, 且覆盖的物种数量(本基因树内)达2/3。
    '''

    t = Tree(each_tree)

    def get_species_in_clade(each_node):
        '''
        获取基因树一个节点后的所有样品名称
        '''
        species_in_clade_list = []
        for each in each_node.get_leaf_names():
            if each.split("++")[0] not in species_in_clade_list:
                species_in_clade_list.append(each.split("++")[0])
        return species_in_clade_list


    def judge1(each_node):
        '''
        输入一个节点, 判断外类群属于该节点后的哪个分支。
        '''
        a = 0
        b = 0
        for each_species in outgroups_list:
            if each_species in str(get_species_in_clade(each_node.children[0])):
                a = 1
            elif each_species in str(get_species_in_clade(each_node.children[1])):
                b = 1
            else:
                pass
        if a == 1 and b == 1: #两个小分支都含有外类群，不符要求
            return None
        elif a == 1 and b == 0: #仅分支1含有外类群
            for each_sample in get_species_in_clade(each_node.children[0]):
                if each_sample.split("++")[0] not in outgroups_list: # 分支1含有不属于外类群的样品
                    return None 
            return 0
        elif a == 0 and b == 1: #仅分支2还有外类群
            for each_sample in get_species_in_clade(each_node.children[1]):
                if each_sample.split("++")[0] not in outgroups_list: # 分支2含有不属于外类群的样品
                    return None
            return 1
        else: #两个分支都不含有外类群
            return None

    def judge2(each_node):
        '''
        输入一个节点, 判断:
        该节点后的所有样品其物种都在ingroups_list中,
        这些物种覆盖了这棵基因树所含的大部分物种(2/3).
        '''
        ingroup_samples_list = []
        for each_sample in get_species_in_clade(each_node):
            if each_sample.split("++")[0] not in ingroups_list:
                return None #有样品不在ingroups_list中
            else:
                if each_sample.split("++")[0] not in ingroup_samples_list:
                    ingroup_samples_list.append(each_sample.split("++")[0])
        all_samples_list = []
        for each_sample in get_species_in_clade(t):
            if each_sample.split("++")[0] in ingroups_list:
                if each_sample.split("++")[0] not in all_samples_list:
                    all_samples_list.append(each_sample.split("++")[0])
        if len(ingroup_samples_list)/len(all_samples_list) > 2/3:
            return True

    calibra_list = []
    n = 0
    for each_node in t.traverse():
        if each_node.children: #如果该节点不是叶子节点
            #判断外类群在该节点后的哪个分支中，如果都没有则遍历下一个节点。
            judge_1 = judge1(each_node)
            if judge_1 != None: #其中一个分支含有外类群
                if judge_1 == 0: #第1个分支中含有外类群
                    if judge2(each_node.children[1]): #第2个分支仅含有内类群，且含有很多
                        calibra_list.append("mrca = a" + str(n) + " " + str(each_node.get_leaf_names()).replace("\'", "").replace(",", "")[1:-1] + "\n")
                        calibra_list.append("min = a" + str(n)  + " " + str(calibrated_age) + "\n")
                        calibra_list.append("max = a" + str(n)  + " " + str(calibrated_age + 0.01) + "\n")
                        n = n + 1
                elif judge_1 == 1: #第2个分支中含有外类群
                    if judge2(each_node.children[0]): #第1个分支仅含有内类群，且含有很多
                        calibra_list.append("mrca = a" + str(n) + " " + str(each_node.get_leaf_names()).replace("\'", "").replace(",", "")[1:-1] + "\n")
                        calibra_list.append("min = a" + str(n)  + " " + str(calibrated_age) + "\n")
                        calibra_list.append("max = a" + str(n)  + " " + str(calibrated_age + 0.01) + "\n")
                        n = n + 1
                else:
                    pass
            else:
                pass
        else:
            pass    
    return calibra_list                   


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
                            "[smooth = 0.0000001]\n"        
                            "randomcv]\n"
                            "cviter = 5\n"
                            "cvsimaniter = 10000\n"
                            "cvstart = 1000\n"
                            "cvstop = 0.0000001\n"
                            "cvmultstep = 0.1\n"
                            "cvoutfile = " + each_tree + ".cvout\n"
                            "outfile = " + each_tree + ".newtree\n")      


for each_tree in get_file_list():
    print(each_tree)
    generate_treePL_conf(each_tree)