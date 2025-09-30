import os
#函数 get_file_list()，得到当前目录中所有的文件名并存入相关列表中
def get_file_list():
  global now_dir
  now_dir = os.getcwd()           #得到当前目录   
  file_temp = os.listdir(path=now_dir)    #得到当前目录中的所有文件名
  global file_name
  file_name = []         #储存基因文件名称
  for each in file_temp:        #将后缀为.fasta的文件加入列表file_name中
      if "RAxML" in each:
          file_name.append(each)

get_file_list()             #得到当前目录中所有的文件名并存入相关列表中
with open(now_dir + "/result.trees", "a") as write_file:
  for each_num in range(100000):
    for each_file in file_name:
      if each_file.split("ortho")[1].split("_cds")[0] == str(each_num):
        with open(now_dir + "/" + each_file, "r") as read_file:
          for each_line in read_file:
            write_file.write(each_line)
          
    









	
           
