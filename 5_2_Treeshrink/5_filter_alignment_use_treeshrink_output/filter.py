import os
from Bio import SeqIO
now_dir = os.getcwd()  
import sys

treethrink_output = "output.txt"

def get_alignment_list():
    # 函数将该文件夹中所有的pep文件加入到一个列表中
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if "cds_maffted.fasta" in each:
            file_list.append(int(each[5:-18]))
    file_list.sort()
    return(file_list)


def main():
    alignment_num_list = get_alignment_list()
    with open(treethrink_output, "r") as read_file:
        n = 0
        for each_line in read_file:
            with open("ortho" + str(alignment_num_list[n]) + "_cds_maffted.fas", "a") as write_file:
                for each_seq in SeqIO.parse("ortho" + str(alignment_num_list[n]) + "_cds_maffted.fasta", "fasta"):
                    if each_seq.id not in each_line:
                        write_file.write(">" + str(each_seq.id) + "\n")
                        write_file.write(str(each_seq.seq) + "\n")
       
            n = n + 1


            
    
main()