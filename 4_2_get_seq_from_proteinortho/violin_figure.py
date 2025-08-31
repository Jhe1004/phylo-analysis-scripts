import os          
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

proteinortho_output_tsv = "myproject.proteinortho.tsv"
species_classify = "violin_species_classify.txt"
max_allow_missing_species_proportion = 0.30
max_allow_low_copy_gene_num = 1

def parsing_proteinortho_output_tsv():
    # 筛选符合最低物种数目的基因出来
    tsv_list = []
    with open(proteinortho_output_tsv, "r") as read_file:
        for each_line in read_file:
            tsv_list.append(each_line)
    species_num = tsv_list[1].count("\t") - 2
    print("species_num:=" + str(species_num))
    min_species_num = species_num - (species_num * max_allow_missing_species_proportion)
    print("min_species_num:=" + str(min_species_num))
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

    # 接下来先将表格转置
    # 先将表格暂时输出成一个tsv文件
    with open("temp_tsv1.tsv", "a") as write_file:
        for each_line in tsv_list:
            write_file.write(each_line)
    # 读取表格，转置，再输出
    tsv = pd.read_csv("temp_tsv1.tsv", sep="\t")
    tsv = tsv.T
    tsv.to_csv("temp_tsv2.tsv", sep="\t")
    #os.remove("temp_tsv1.tsv")

    # 物种名改为分类名称并统计基因数目
    # 先读取物种分类别的txt
    species_classify_list = []
    with open(species_classify, "r") as read_file:
        for each_line in read_file:
            if len(each_line) > 2:
                if each_line[-1] == "\n":
                    species_classify_list.append(each_line[:-1].split("\t"))
    # 按照species_classify_list中的分类修改列表
    with open("temp_tsv3.tsv", "a") as write_file:
        write_file.write("sp\tX\tY\n")
        with open("temp_tsv2.tsv", "r") as read_file:
            for each_line in read_file:
                gene_number = each_line.count("\t") - each_line.count("*")
                if each_line.split(".fasta")[0] in species_classify_list[0]:
                    write_file.write(each_line.split(".fasta")[0] + "\tLN\t" + str(gene_number) + "\n")
                elif each_line.split(".fasta")[0] in species_classify_list[1]:
                    write_file.write(each_line.split(".fasta")[0] + "\tSG\t" + str(gene_number) + "\n")
                else:
                    pass
    #os.remove("temp_tsv2.tsv")

def violin():
    violin_data = pd.read_csv("temp_tsv3.tsv", sep="\t")
    f, ax = plt.subplots(figsize=(20, 10))
    sns.violinplot(x="X",y="Y",data=violin_data)
    sns.stripplot(x="X", y="Y", data=violin_data, color="black", alpha=.6, jitter=0.02)
    plt.savefig("violin.png", dpi = 120)
    





    


parsing_proteinortho_output_tsv()    
violin()
#os.remove("temp_tsv3.tsv")