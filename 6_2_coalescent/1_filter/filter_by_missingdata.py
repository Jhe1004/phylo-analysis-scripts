'''
并联法在astral之前的筛选，通过之前raxml输出文件中的missing data数量进行筛选
输入文件有alignment排序文件，以及RAxML_bipartitions.*树文件，这两种文件需要放置在本文件夹中
需要设置允许的最大missing data阈值以及输出文件夹的名称，输出文件夹中将包含筛选得到的alignment排序文件和树文件
'''
import os
from Bio import SeqIO
now_dir = os.getcwd() 

max_missingdata = 0.2   #设置允许的最小bootstrap值
output_dir = "max_missingdata_02"

def get_file_list():
    #获取当前文件夹中指定文件
    #输入：无
    #输出：指定文件列表：file_list
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if ".fasta" in each:
            file_list.append(each)
    return file_list

def main():
    #函数mian: 主控函数
    #输入：无
    #输出：脚本运行结果
    os.mkdir(output_dir)
    file_list = get_file_list()
    file_dict = {} 
    species_number_max = 0
    for each_file in file_list:
        file_dict[each_file] = SeqIO.to_dict(SeqIO.parse(each_file, "fasta"))
        if len(file_dict[each_file]) >= species_number_max:
            species_number_max = len(file_dict[each_file])
    for each_file in file_dict:
        if len(file_dict[each_file])/species_number_max >= 1 - max_missingdata:
            os.system("cp " + each_file + " ./" + output_dir)
            os.system("cp RAxML_bipartitions." + each_file[:-6] + " ./" + output_dir)
            


  


   
main()