'''
计算两个物种之间的ks值
'''
import os
import re
from Bio import SeqIO
now_dir = os.getcwd()    

species1 = "Hydrastis_canadensis_2.fasta.transdecoder.pep"
species2 = "Anemone_hupehensis_18.fasta.transdecoder.pep"

#函数get_file_list：获取当前文件夹中指定文件
def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if ".fasta" in each:
            file_list.append(each)
    return file_list

def codeml():
    with open("codeml.ctl", "a") as write_file:
        write_file.write("seqfile = " + "temp.phy" + "\n")
        write_file.write("treefile = stewart.trees" + "\n")
        write_file.write("outfile = outfile" + "\n")
        write_file.write("noisy = 0" + "\n")
        write_file.write("verbose = 0" + "\n")
        write_file.write("runmode = -2" + "\n")
        write_file.write("seqtype = 1" + "\n")
        write_file.write("CodonFreq = 2" + "\n")
        write_file.write("clock = 0" + "\n")
        write_file.write("aaDist = 0" + "\n")
        write_file.write("model = 0" + "\n")
        write_file.write("NSsites = 0  " + "\n")
        write_file.write("icode = 0" + "\n")
        write_file.write("Mgene = 0" + "\n")
        write_file.write("fix_kappa = 0" + "\n")
        write_file.write("kappa = 2" + "\n")
        write_file.write("fix_omega = 0" + "\n")
        write_file.write("omega = .4" + "\n")
        write_file.write("fix_alpha = 1" + "\n")
        write_file.write("alpha = 0" + "\n")
        write_file.write("Malpha = 0" + "\n")
        write_file.write("ncatG = 8" + "\n")
        write_file.write("getSE = 0" + "\n")
        write_file.write("RateAncestor = 1" + "\n")
        write_file.write("Small_Diff = .5e-6" + "\n")
        write_file.write("cleandata = 1" + "\n")
        write_file.write("method = 0" + "\n")
    os.system("codeml")
    with open("result.txt","a") as write_file:
        with open("2ML.dS","r") as read_file:
            read_file.readline()
            read_file.readline()
            ks = read_file.readline()
            print(ks.split(" ")[-1])
            write_file.write(ks.split(" ")[-1])
    os.remove("2ML.dN")
    os.remove("2ML.dS")
    os.remove("2ML.t")
    os.remove("2NG.dN")
    os.remove("2NG.dS")
    os.remove("2NG.t")
    os.remove("codeml.ctl")
    os.remove("outfile")
    os.remove("rst")
    os.remove("rst1")
    os.remove("rub")    

#将直系同源基因文件中需要计算的基因对提取出来
def calculate_ks(each_fasta):
    file_dict = {}
    file_dict[each_fasta] = SeqIO.to_dict(SeqIO.parse(each_fasta, "fasta"))
    if species1 in file_dict[each_fasta] and species2 in file_dict[each_fasta]:
        with open("temp.phy", "a") as write_file:
            write_file.write(" 2 " + str(len(file_dict[each_fasta][species1].seq)) + "\n")
            write_file.write(str(file_dict[each_fasta][species1].id) + "  " + str(file_dict[each_fasta][species1].seq) + "\n")
            write_file.write(str(file_dict[each_fasta][species2].id) + "  " + str(file_dict[each_fasta][species2].seq) + "\n")
        codeml()
        os.remove("temp.phy")
  

def main():
    fasta_file_list = get_file_list()
    for each_fasta in fasta_file_list:
        calculate_ks(each_fasta)

if __name__ == "__main__":
    main()
