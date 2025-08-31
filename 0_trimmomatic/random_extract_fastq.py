import os           #导入os模块
import random
now_dir = os.getcwd()           #
file_name = os.listdir(path=now_dir)
file_name2 = []



def build_dict():
    with open(now_dir + "/" + fastq_name[:-2], "w") as write_file:
        with open(now_dir + "/" + fastq_name, "r") as read_file:
            for each_line in read_file:
                if each_line:   
                    if each_line[0] == "@":
                        random_num = random.random()
                        if random_num < 0.05:
                            write_file.write(each_line)
                            a = read_file.readline()
                            write_file.write(a)
                            b = read_file.readline()
                            write_file.write(b)
                            c = read_file.readline()
                            write_file.write(c)
            
    
        
for each in file_name:
    if "SRR7271011_1" in each:
        file_name2.append(each)

for fastq_name in file_name2:
    build_dict()

