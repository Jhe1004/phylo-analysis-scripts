import os
import re
from Bio import SeqIO


def get_file_list():
  global now_dir
  now_dir = os.getcwd()          
  file_temp = os.listdir(now_dir)    
  global file_name
  file_name = []       
  for each in file_temp:        
      if ".txt" in each:
          file_name.append(each)


def Reduce_01_5(fasta_name):
  with open(now_dir + "/" + fasta_name + "new", "a") as write_file:
    with open(now_dir + "/" + fasta_name, "r") as read_file:
      for each in read_file:
        try:
          if float(each[:-1]) >= 0.1:
            if float(each[:-1]) <= 5:
              write_file.write(each)
        except:
          pass
        

get_file_list()          

for fasta_name in file_name:
  Reduce_01_5(fasta_name)

