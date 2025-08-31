from multiprocessing import Pool
import os
import subprocess
from ete3 import Tree

def get_file_list():      
    file_name = []       
    for each in os.listdir(os.getcwd()):        
        if ".fasta" in each:
            file_name.append(each)
    return file_name



def get_ks_from_raxml(each_fasta):
    #subprocess.call("raxmlHPC-PTHREADS -T 1 -n "+ each_fasta.replace(".fasta", "") + " -s " + each_fasta.replace(".fasta", ".fas") + " -m GTRGAMMA -p 12345 -f a -N 2 -x 12345", shell = True)
    besttree = Tree("RAxML_bestTree." + each_fasta.replace(".fasta", ""))
    print(besttree.get_ascii(show_internal=True))
    internal_node_name = n
    for node in besttree.traverse():
        if node.is_leaf():
            pass
        else:
            node.add_features(name = str(internal_node_name))
            internal_node_name = internal_node_name + 1
    print(besttree.get_ascii(show_internal=True))
    




file_list = get_file_list()
for each_fasta in file_list:
    with open(each_fasta, "r") as read_file:
        n = 0
        for each_line in read_file:
            if each_line[0] == ">":
                n = n + 1
    if n < 4:
        pass
    else:
        get_ks_from_raxml(each_fasta)
            
