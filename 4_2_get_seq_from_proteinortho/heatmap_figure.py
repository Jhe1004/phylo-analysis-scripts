import os          
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

'''
python get_seq_from_proteinortho.py -i myproject.proteinortho.tsv -o output2 -p 0.3 -c 1 -n name_list.txt
'''

#解析参数
#参数分为必须参数(required)和可选参数(additional)
parser = argparse.ArgumentParser(description="Options for get_seq_from_proteinortho.py", 
                                add_help=True)
required = parser.add_argument_group("Required arguments")
required.add_argument('-i', '--infile', action="store", metavar='\b', 
                      type=str, required=True, default="myproject.proteinortho.tsv", 
                      help="Name of the hyde proteinortho file")  
required.add_argument('-p', '--missing', action="store", metavar='\b', type=float, 
                      required=True, default=0.3, help="max allow missing species proportion")  
required.add_argument('-c', '--copy', action="store", metavar='\b', type=int, 
                      required=True, default=1, help="max allow low copy num")
additional = parser.add_argument_group("Additional arguments")
additional.add_argument('-n', '--name_order', action="store", metavar='\b', type=str, 
                        help="The order in which species names appear on a heat map")                                    

args                                        = parser.parse_args()
proteinortho_output_tsv                     = args.infile
max_allow_missing_species_proportion        = args.missing
max_allow_low_copy_gene_num                 = args.copy
species_order                               = args.name_order



def parsing_proteinortho_output_tsv():
    # 筛选符合最低物种数目的基因出来
    tsv_list = []
    with open(proteinortho_output_tsv, "r") as read_file:
        for each_line in read_file:
            tsv_list.append(each_line)
    species_num = tsv_list[1].count("\t") - 2
    min_species_num = species_num - (species_num * max_allow_missing_species_proportion)
    tsv_list_temp = []
    tsv_list_temp.append("\t".join(tsv_list[0].split("\t")[3:]))
    for each in tsv_list[1:]:
        if int(each.split("\t")[0]) >= min_species_num:
            each = "\t".join(each.split("\t")[3:])
            tsv_list_temp.append(each)
    tsv_list = tsv_list_temp

    # 筛选不超过低拷贝基因数目的基因出来
    tsv_list_temp = []
    tsv_list_temp.append(tsv_list[0])
    for each in tsv_list[1:]:
        max_gene_number = 0
        for each_gene in each.split("\t"):
            gene_copy_number = each_gene.count(",") + 1
            if gene_copy_number > max_allow_low_copy_gene_num:
                max_gene_number = gene_copy_number
        if max_gene_number > max_allow_low_copy_gene_num:
            pass
        else:
            tsv_list_temp.append(each)
    tsv_list = tsv_list_temp

    # 转化为“0”“1”矩阵，没有基因为“0”，有基因为“1”
    tsv_list_temp = []
    tsv_list_temp.append(tsv_list[0])
    for each in tsv_list[1:]:
        each_tem = []
        for each_gene in each.split("\t"):
            if each_gene == "*" or each_gene == "*\n":
                each_tem.append("0")
            else:
                each_tem.append("1")
        tsv_list_temp.append("\t".join(each_tem) + "\n")
    tsv_list = tsv_list_temp

    # 接下来先将表格转置，之后让物种按照一个特殊的给定的排列方式排列
    # 先将表格暂时输出成一个tsv文件
    with open("temp_tsv1.tsv", "a") as write_file:
        for each_line in tsv_list:
            write_file.write(each_line)
    # 读取表格，转置，再输出
    tsv = pd.read_csv("temp_tsv1.tsv", sep="\t")
    tsv = tsv.T
    tsv.to_csv("temp_tsv2.tsv", sep="\t")
    # 再读取，排列物种名，再输出
    if species_order:
        tsv_list = []
        with open("temp_tsv2.tsv", "r") as read_file:
            for each_line in read_file:
                tsv_list.append(each_line)
        tsv_list_temp = []
        with open(species_order, "r") as read_file:
            for each_line in read_file:
                if len(each_line) > 2:
                    have_species = 0
                    for each_specie in tsv_list:
                        if each_line.replace("\n", "") == each_specie.split("\t")[0]:
                            tsv_list_temp.append(each_specie)
                            have_species = 1
                    if have_species == 0:
                        print("the species name " + each_line[:-1] + " not in the proteinortho output")
        with open("temp_tsv.tsv", "a") as write_file:
            for each_line in tsv_list_temp:
                write_file.write(each_line)

    os.remove("temp_tsv1.tsv")
    #os.remove("temp_tsv2.tsv")
      

def heatmap():
    if species_order:
        tsv = pd.read_csv("temp_tsv.tsv", sep="\t", header=None, index_col=0)
        os.remove("temp_tsv.tsv")
    else:
        tsv = pd.read_csv("temp_tsv2.tsv", sep="\t", index_col=0)
    f, ax = plt.subplots(figsize=(60, 40))
    sns.heatmap(tsv, ax=ax, cmap='GnBu')
    plt.savefig("heatmap.png", dpi = 120)

parsing_proteinortho_output_tsv()
heatmap()
os.remove("temp_tsv2.tsv")