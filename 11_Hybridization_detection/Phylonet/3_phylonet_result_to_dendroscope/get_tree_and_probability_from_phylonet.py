import dendropy
import os
import re
now_dir = os.getcwd()           #得到当前目录  


def get_file_list():
  #获取当前文件夹中指定文件
  #输入：无
  #输出：指定文件列表：file_list
  file_temp = os.listdir()
  file_list = []
  for each in file_temp:       
      if ".txt" in each:
          file_list.append(each)
  return file_list   


def get_tree(each_file):
    with open(each_file, "r") as read_file:
        for each_line in read_file:
            if "Visualize in Dendroscope" in each_line:
                with open(now_dir + os.sep + "tree" + os.sep + each_file, "a") as write_file:
                    write_file.write(each_line[27:])
                
def get_probability(each_file):
    with open(now_dir + os.sep + "probability" + os.sep + "result.txt", "a") as write_file:
        with open(each_file, "r") as read_file:
            for each_line in read_file:
                if "Total log probability" in each_line:                 
                    write_file.write(each_file + "," + each_line[23:])
                    break




def main(file_list):
    for each_file in file_list:
        get_tree(each_file)
        get_probability(each_file)





file_list = get_file_list()
main(file_list)