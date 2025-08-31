import os          
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

proteinortho_output_file = "myproject.proteinortho.tsv"
output_seq_dir = "output"


'''
首先从proteinortho的结果中提取所有基因家族文件
'''
# make output dir
try:
    os.makedirs(os.getcwd() + os.sep + output_seq_dir)
except:
    print("The output directory \"" + output_seq_dir + "\" already exists, script stop.")
    sys.exit()

# Import the FASTA sequence into memory
file_dict = {}  
for each in os.listdir(os.getcwd()):        
    if ".pep" in each:
        file_dict[each] = SeqIO.to_dict(SeqIO.parse(each, "fasta"))
    elif ".cds" in each:
        file_dict[each] = SeqIO.to_dict(SeqIO.parse(each, "fasta"))
    else:
        pass

# get sequences from ".cds" or ".pep" files.            
def get_nucleotide(IDs): 
    return str(file_dict[IDs.split("++")[0][:-3] + "cds"][IDs.split("++")[1]].seq)     
def get_protein(IDs): 
    return str(file_dict[IDs.split("++")[0][:-3] + "pep"][IDs.split("++")[1]].seq)    

# main function
def get_seq():
    protein_table = pd.read_csv(proteinortho_output_file, sep="\t")
    species_num = protein_table.shape[1] - 3
    columns_values = protein_table.columns.values
    i = 0
    for index, row in protein_table.iterrows():
        i = i + 1           
        if row.iloc[0] == row.iloc[1]:
            pass
        elif row.iloc[1] < 4:
            pass
        else:
            nucleotide_fasta = open(os.getcwd() + os.sep + output_seq_dir + os.sep + "ortho" + str(i) + "_cds.fasta", "a")
            protein_fasta = open(os.getcwd() + os.sep + output_seq_dir + os.sep + "ortho" + str(i) + "_pep.fasta", "a")
            for n in range(3, protein_table.shape[1]):
                if row.iloc[n] != "*":
                    for each_seq in row.iloc[n].split(","):
                        nucl_seq = get_nucleotide(columns_values[n] + "++" + each_seq)
                        nucleotide_fasta.write(">" + columns_values[n] + "++" + each_seq + "\n" + nucl_seq + "\n")
                        prot_seq = get_protein(columns_values[n] + "++" + each_seq)
                        protein_fasta.write(">" + columns_values[n] + "++" + each_seq + "\n" + prot_seq + "\n")       
                else:
                    pass
            
get_seq()