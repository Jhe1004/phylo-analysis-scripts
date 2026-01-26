import os
import re
from Bio import SeqIO
now_dir = os.getcwd()    

def fasta2dict(fasta_file):
    #将fasta文件转化成一个python列表
    res_dict = {}
    with open(fasta_file, "r") as read_file:
        for each_line in read_file:
            if each_line[0] == ">":
                seq_name = each_line.split(" ")[0][1:].replace("\n", "").replace("/","_").replace("\\","_")
                res_dict[seq_name] = ""
            else:
                res_dict[seq_name] = res_dict[seq_name] + each_line.replace("\n", "")
    return res_dict, len(res_dict[seq_name])
    


def get_file_list():
    #函数get_file_list：获取当前文件夹中指定文件
    #输入：无
    #输出：指定文件列表：file_list
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:
        if ".fas" in each:
            file_list.append(each)
    return file_list


def main():
    fasta_file_list = get_file_list()
    concat_list = []
    name_list = []
    for each_file in fasta_file_list:
        fasta_dict = fasta2dict(each_file)
        concat_list.append(fasta_dict)
        for each_name in fasta_dict[0]:
            if each_name not in name_list:
                name_list.append(each_name)
    for each_name in name_list:
        with open(each_name + ".temp", "a") as write_file:
            write_file.write(">" + each_name + "\n")
            for each_dict in concat_list:
                if each_name in each_dict[0]:
                    write_file.write(each_dict[0][each_name])
                else:
                    write_file.write("?"*each_dict[1])
            write_file.write("\n")
    os.system("cat *.temp > result.fasta")
    os.system("rm *.temp")
            



if __name__ == "__main__":
    main()
    