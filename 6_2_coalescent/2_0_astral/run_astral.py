'''
直接运行astral，
输入文件是raxml输出的树文件,且合并到一起
'''
import os
now_dir = os.getcwd() 

input_trees_dir = "input"

def get_file_list():
    #获取当前文件夹中指定文件
    #输入：无
    #输出：指定文件列表：file_list
    file_temp = os.listdir(now_dir + os.sep + input_trees_dir)
    file_list = []
    for each in file_temp:       
        if ".trees" in each:
            file_list.append(each)
    return file_list

def main():
    trees_list = get_file_list()
    for each_trees in trees_list:
        os.system("java -jar astral.5.7.8.jar -t 2 --outgroup Symphalangus.fasta.transdecoder.pep -i " + input_trees_dir + "/" + each_trees + " -o " + each_trees.replace("tree","") + ".tree")


main()



  
		
           
