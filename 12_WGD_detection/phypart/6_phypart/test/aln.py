import os
import re
now_dir = os.getcwd()    

#函数get_file_list：获取当前文件夹中指定文件
def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if ".tree" in each:
            file_list.append(each)
    return file_list




def main():
    fasta_file_list = get_file_list()
    for each_tree_file in fasta_file_list:
        with open(each_tree_file, "r") as read_file:
            for each_line in read_file:
                #new_line = re.sub("\+\+.*?\.p[0-9]", "", each_line, count=0, flags=0)
                new_line = re.sub("\.fasta\.transdecoder\.pep\+\+.*?\.p[0-9]", "", each_line, count=0, flags=0)
            with open(each_tree_file[:-10] + ".tree", "a") as write_file:
                write_file.write(new_line)

if __name__ == "__main__":
    main()



