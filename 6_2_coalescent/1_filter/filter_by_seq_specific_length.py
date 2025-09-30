'''
并联法在astral之前的筛选，通过之前raxml建树时的alignment长度进行筛选
输入文件有alignment排序文件，以及RAxML_bipartitions.*树文件，这两种文件需要放置在本文件夹中
需要设置允许的最小长度阈值以及输出文件夹的名称，输出文件夹中将包含筛选得到的alignment排序文件和树文件
'''
import os
from Bio import SeqIO
now_dir = os.getcwd()

min_seq_len = 500   #设置允许的序列最小长度
max_seq_len = 1000
output_dir = "seq_length_500_1000"

def get_file_list():
  #函数get_file_list：获取当前文件夹中指定文件
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
  fasta_file_list = get_file_list()
  for each_file in fasta_file_list:
    for each_seq in SeqIO.parse(each_file, "fasta"):
      seq_len = len(str(each_seq.seq))
      break
    if seq_len > min_seq_len and seq_len < max_seq_len:
      os.system("cp " + each_file + " ./" + output_dir)
      #os.system("cp RAxML_bipartitions." + each_file[:-6] + " ./" + output_dir)


main()