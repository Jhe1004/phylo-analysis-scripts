import os
import re

'''
函数 get_file_list: 获取当前文件夹中符合目标扩展名的文件
输入 无，将本脚本放置在目标文件夹中即可
输出 file_name：所有文件的名称列表
'''
def get_file_list():      
    file_name = []       
    for each in os.listdir(os.getcwd()):        
        if ".fas" in each:
            file_name.append(each)
    return file_name

file_name = get_file_list()
for each_name in file_name:
    os.rename(each_name, each_name[:-16] + ".fasta")

'''
file_name = get_file_list()
for each in file_name:
    with open(each[:-6] + ".fas", "a") as write_file:
        with open(each, "r") as read_file:
            for each_line in read_file:
                if each_line[0] == ">":

                    pattern1 = ">ortho[0-9]*\|"
                    newseq = re.sub(pattern1, ">", each_line)

                    pattern1 = "\.fas.*"
                    newseq = re.sub(pattern1, "", newseq)
                  

                    write_file.write(newseq)                    
                else:
                    write_file.write(each_line)


'''



