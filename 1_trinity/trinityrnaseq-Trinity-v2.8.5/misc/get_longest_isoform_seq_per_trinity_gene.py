import os
def get_file_list():
  global now_dir
  now_dir = os.getcwd()          
  file_temp = os.listdir(now_dir)    
  global file_name
  file_name = []       
  for each in file_temp:        
      if ".fasta" in each:
          file_name.append(each)






get_file_list()          

for fasta_name in file_name:   
    print(fasta_name)
    os.system("perl get_longest_isoform_seq_per_trinity_gene.pl "+ fasta_name + " > " + fasta_name[:-2]) 
           
