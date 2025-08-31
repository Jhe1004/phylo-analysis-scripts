import os
from ete3 import Tree
now_dir = os.getcwd()    

outgroups = "Anemoclema_glauciifolium_19.fasta.transdecoder.pep"


#函数get_file_list：获取当前文件夹中指定文件
def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if "RAxML" in each:
            file_list.append(each)
    return file_list


def reroot(each_trees_file):
    with open(each_trees_file + "_reroot.trees", "a") as write_file:
        with open(each_trees_file, "r") as read_file:
            for each_line in read_file:
                if len(each_line) > 2:
                    t = Tree(each_line)
                    t.set_outgroup(outgroups)
                    write_file.write(t.write() + "\n")


def main():
    fasta_file_list = get_file_list()
    for each_trees_file in fasta_file_list:
        reroot(each_trees_file)

if __name__ == "__main__":
    main()