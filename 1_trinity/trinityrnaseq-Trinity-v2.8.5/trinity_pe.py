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


for each in file_name:
    command = "Trinity --seqType fq --min_kmer_cov 2 --left " + each + "_1.fq --right " + each + "_2.fq  --output " + each + "_trinity --max_memory 60G --CPU 23 --full_cleanup"
    print("run command: " + command)
    os.system(command)         
