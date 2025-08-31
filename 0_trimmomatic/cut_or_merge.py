import os   
del_gene_list = []

def get_file_list():
  global now_dir
  now_dir = os.getcwd()         
  file_temp = os.listdir(path=now_dir)   
  global file_name
  file_name = []   
  for each in file_temp:    
      if ".fq" in each:
          file_name.append(each)
          
          
          
def del_gene():
    with open(now_dir + "/" + fasta_name[:-2] + "fastq", "a") as write_file:
        with open(now_dir + "/" + fasta_name, "r") as read_file:
            n = 1
            for each_line in read_file:
                if each_line[0] == "@":
                    n  = n + 1
                    if n % 6 == 0:
                        write_file.write(each_line)
                        write_file.write(read_file.readline())
                        write_file.write(read_file.readline())
                        write_file.write(read_file.readline())
                    
                
    
          
            
          
def merge_gene(fasta_name):
    with open(now_dir + "/result.fastq", "a") as write_file:
        with open(now_dir + "/" + fasta_name, "r") as read_file:
            for each_line in read_file:
                write_file.write(each_line)          
          
          
          
          
          
          
          
          
          
          
          
          
          
          

get_file_list()     


for fasta_name in file_name:
    merge_gene(fasta_name)  
    
