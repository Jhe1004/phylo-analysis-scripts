'''
从raxml的结果中筛选出phylonetworks需要的输入树，长度太多或者没有包含所有要求的物种树的
输入：raxml当时输入的alignment文件，用于判断建树时矩阵的长度
raxml输出的bootstrap文件
一个包含带提取物种名称的文件，文件中每行代表一个物种
输出：筛选好的bootstrap文件
bslistfiles文件，文件中每行是一个bootstrap文件名
'''

import os
from Bio import SeqIO
from ete3 import Tree
now_dir = os.getcwd()    

need_species_list_file = "result.txt"  #该文件每一行写入一个待提取的物种名称
Sequence_minimum_length = 1000 #设置允许的最小建树矩阵长度


#函数get_file_list：获取当前文件夹中指定文件
def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if ".fasta" in each:
            file_list.append(each)
    return file_list


def long_enough(each_fasta):
    for each_seq in SeqIO.parse(each_fasta, "fasta"):
        if len(str(each_seq.seq)) > 1000:
            return True
        else:
            return False
        break

def write_bslist():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if ".tree" in each:
            file_list.append(each)
    for each in file_list:
        with open("BSlistfiles", "a") as write_file:
            write_file.write(each + "\n")

def main():
    fasta_file_list = get_file_list()
    species_list = []
    with open(need_species_list_file, "r") as read_file:
        for each_line in read_file:
            if len(each_line) > 2:
                species_list.append(each_line.replace("\n",""))
    for each_fasta in fasta_file_list:
        if long_enough(each_fasta):
            write_file = open("RAxML_bootstrap." + each_fasta.replace(".fasta", "") + "_new.trees", "a")
            with open("RAxML_bootstrap." + each_fasta.replace(".fasta", ""), "r") as read_file:
                for each_line in read_file:
                    if len(each_line) > 2:
                        t = Tree(each_line.replace("\n",""))
                        try:
                            t.prune(species_list)
                            write_file.write(t.write() + "\n")
                        except:
                            write_file.close()
                            os.remove("RAxML_bootstrap." + each_fasta.replace(".fasta", "") + "_new.trees")
                            break
            write_file.close()
    write_bslist()
                        

                        
            


if __name__ == "__main__":
    main()
