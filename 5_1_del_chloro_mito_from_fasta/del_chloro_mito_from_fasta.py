import os
import re
from Bio import SeqIO
now_dir = os.getcwd()    


#函数get_gb_file：文件夹中的gb文件
def get_gb_file():
    file_temp = os.listdir(now_dir)
    gb_list = []
    for each in file_temp:       
        if each[-2:] == "gb":
            gb_list.append(now_dir + os.sep + each)
    return gb_list

#函数get_pep_file：获取pep文件
def get_pep_file():
    file_temp = os.listdir(now_dir)
    pep_list = []
    for each in file_temp:       
        if each[-9:] == "pep.fasta":
            pep_list.append(now_dir + os.sep + each)
    return pep_list  

#函数get_cds_from_ref_gb：将gb文件中的所有编码区全部提取出来
def get_protein_seq_from_ref_gb(each_gb):
    with open("ref.protein", "a") as write_file:
        for record in SeqIO.parse(each_gb, "genbank"):
            tmp_list = [] 
        for gene in record.features:
            if gene.type == "CDS":
                if str(gene.qualifiers["gene"]) not in tmp_list:
                    tmp_list.append(str(gene.qualifiers["gene"]))
                    write_file.write(">" + str(gene.qualifiers["gene"])[2:-2] + "\n")
                    write_file.write(str(gene.qualifiers["translation"])[2:-2] + "\n")
            else:
                pass

def get_fasta_seq_number(each_pep):
    n = 0
    with open(each_pep, "r") as read_file:
        for each_line in read_file:
            if each_line[0] == ">":
                n = n + 1
    return n




def main():
    for each_gb in get_gb_file():
        get_protein_seq_from_ref_gb(each_gb)
    os.system("makeblastdb -in ref.protein -input_type fasta -title ref -out ref -dbtype prot")
    with open("del_file.txt", "a") as write_file:
        for each_pep in get_pep_file():
            os.system("blastp -task blastp-short -query " + each_pep + " -db ref -out ref.result -evalue 0.00001 -outfmt \"7 std qlen slen\"")
            fasta_seq_number = get_fasta_seq_number(each_pep)
            with open("ref.result", "r") as read_file:
                n = 0
                for each_line in read_file:
                    if each_line == "# 0 hits found\n":
                        n = n + 1
                if (fasta_seq_number - n) / fasta_seq_number >= 0.1:
                    os.remove(each_pep)
                    os.remove(each_pep.replace("pep", "cds"))
                    write_file.write(each_pep + "\n" + each_pep.replace("pep", "cds") + "\n")
                else:
                    pass
            os.remove("ref.result")



main()
os.system("rm ref.*")