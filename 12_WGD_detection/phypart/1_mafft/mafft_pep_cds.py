'''
该脚本需要两个输入文件，分别是proteinortho输出的pep(蛋白质)文件以及proteinortho同时输出的cds(碱基)文件。脚本会先排序蛋白质文件，之后根据蛋白质文件再排序碱基文件。
脚本会同时输出蛋白质文件和碱基文件。
'''

from multiprocessing import Pool
import os, time, random
from Bio import SeqIO, Seq
import argparse
import sys


#解析参数

parser = argparse.ArgumentParser(description="Options for mafft_pep_cds.py",
                                            add_help=True)
additional = parser.add_argument_group("additional arguments")    
additional.add_argument('-n', '--num_cpu', action="store", type=int, default=60,
                        metavar='\b', help="The maximum number of CPUs that this script can be used (Two usable CPUs means that the\
                             script will analyze two matrices at the same time), default = 1")
additional.add_argument('-t', '--thread', action="store", type=int, default=1,
                        metavar='\b', help="Number of threads for each analysis of matrices, default = 1, (If only a single or a\
                             small number of matrices are analyzed, the setting of the '-n' parameter will not increase the computational efficiency,\
                              then you can choose specify this parameter to increase the computational efficiency. but if you need to align a\
                             large number of short gene files (such like preparing data for ASTRAL analysis), it is highly recommended to\
                                  set this parameter to 1 and set the '-n' to the maximum allowable number of CPUs\
                             for your computer")
additional.add_argument('-j', '--just2mafft', action="store_true",
                        help="just run mafft for each fasta")  

args = parser.parse_args()
num_cpu    = args.num_cpu
thread     = args.thread
just2mafft     = args.just2mafft



#定义需要的进程数th

th = num_cpu



#函数 get_file_list: 获取当前文件夹中符合目标扩展名的文件
#输入 无，将本脚本放置在目标文件夹中即可
#输出 file_name：所有文件的名称列表

def get_file_list():      
    file_name = []       
    for each in os.listdir(os.getcwd()):        
        if "pep.fas" in each:
            file_name.append(each)
    return file_name


#函数 get_th_list: 将所有文件平均的分给各个进程
#输入 file_name: 所有文件的名称列表
#输出 th_list： 一个列表，列表中的各个元素为各个进程所应该处理的文件

def get_th_list(file_name):
    th_list = []
    for each_num in range(th):
        th_list.append([])
    print("本程序共使用 " + str(len(th_list)) + " 个进程")
    n = 0
    for each_file in file_name:       
        th_list[n].append(each_file)
        if n == th - 1:
            n = 0
        else:
            n = n + 1
    return th_list




#函数 aa2nt_aln: mafft蛋白文件之后，根据蛋白质文件来排序碱基文件
#输入 aa_aln: 排序好的蛋白质文件；nt_fasta：未排序的碱基文件；outfile：输出文件名称
#输出 无

def aa2nt_aln(aa_aln,nt_fasta,outfile):
    # store the nucleotide sequences in a dictionary
    nt_file = open(nt_fasta,"r")
    nt_dict = SeqIO.to_dict(SeqIO.parse(nt_file, "fasta"))
    nt_file.close()
    # store the nucleotide sequences in a dictionary
    aa_file = open(aa_aln,"r")
    aa_dict = SeqIO.to_dict(SeqIO.parse(aa_file, "fasta"))
    aa_file.close()
    outfile = open(outfile,"a")
    # read through an aa sequence one site at a time
    # if the site is not a gap insert the corresponding codon into a new nt sequence
    # if it is a gap, insert three gap characters
    for seq in aa_dict:
        new_seq=""
        counter=0
        for character in aa_dict[seq]:
            if character != '-':
                if character != '*':
                    new_seq = new_seq+nt_dict[seq].seq[counter:counter+3]
                    counter = counter+3
                else:
                    new_seq = new_seq+"---"
                    counter = counter+3					
            else:
                new_seq = new_seq+"---"
        outfile.write(">" + str(seq) + "\n")
        outfile.write(str(new_seq) + "\n") 




def main_software(each_file_name_list):
    print(each_file_name_list)
    for each in each_file_name_list:
        #print("mafft --thread " + str(thread) + " " + os.getcwd() + os.sep + each + " > " + os.getcwd() + os.sep + each[:-4] + "_maffted.fasta")
        os.system("mafft --thread " + str(thread) + " " + os.getcwd() + os.sep + each + " > " + os.getcwd() + os.sep + each[:-6] + "_maffted.fasta")
        #print(each[:-6] + "_maffted.fasta", each[:-9] + "cds.fasta", each[:-2])
        
        if not just2mafft:
            aa2nt_aln(each[:-6] + "_maffted.fasta", each[:-9] + "cds.fasta", each[:-9] + "cds_maffted.fasta")
        else:
            os.system("mafft --thread " + str(thread) + " " + os.getcwd() + os.sep + each[:-9] + "cds.fasta" + " > " + os.getcwd() + os.sep + each[:-9] + "cds_maffted.fasta")

file_name = get_file_list()
th_list = get_th_list(file_name)

p = Pool(th)
for i in range(th):
    p.apply_async(main_software,(th_list[i],))

print("----start----")
p.close()  # 关闭进程池，关闭后po不再接收新的请求
p.join()  # 等待po中所有子进程执行完成，再执行下面的代码,可以设置超时时间join(timeout=)
print("-----end-----") 
  
		
           
