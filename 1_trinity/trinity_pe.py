import os
#函数 get_file_list()，得到当前目录中所有的文件名并存入相关列表中
def get_file_list():
	global now_dir
	now_dir = os.getcwd()           #得到当前目录   
	file_temp = os.listdir(path=now_dir)    #得到当前目录中的所有文件名
	global file_name
	file_name = []         #储存基因文件名称
	for each in file_temp:        #将后缀为.fasta的文件加入列表file_name中
		if ".fq" in each:
			if each[:-5] not in file_name:
				file_name.append(each[:-5])

get_file_list()             #得到当前目录中所有的文件名并存入相关列表中

now_dir = os.getcwd()
for each in file_name:
    command = "singularity exec -e trinityrnaseq.v2.15.1.simg Trinity --seqType fq --left " + now_dir + os.sep + each + "_1.fq --right " + now_dir + os.sep + each + "_2.fq --max_memory 40G --CPU 40 --output " + now_dir + os.sep + each + "_trinity --full_cleanup"
    print("run command: " + command)
    os.system(command)         
