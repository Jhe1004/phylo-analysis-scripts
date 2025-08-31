'''
直接运行astral，
输入文件是raxml输出的RAxML_bipartitions.*树文件
需要设置输入树所在的文件夹名称(默认input_trees)以及输出树的名称(默认output_tree)
'''
import os
now_dir = os.getcwd() 

input_trees_dir = "input_trees"
output_tree = "result.trees"

def get_file_list():
    #获取当前文件夹中指定文件
    #输入：无
    #输出：指定文件列表：file_list
    file_temp = os.listdir(now_dir + os.sep + input_trees_dir)
    file_list = []
    for each in file_temp:       
        if "RAxML" in each:
            file_list.append(each)
    return file_list

def main():
    tree_list = get_file_list()
    with open(output_tree + "_temp_treelist", "a") as write_file:
        for each_tree in tree_list:
            with open(now_dir + os.sep + input_trees_dir + os.sep + each_tree, "r") as read_file:
                write_file.write(read_file.readline())
    os.system("java -jar astral.jar -i " + output_tree + "_temp_treelist -o " + output_tree + ".tree")


main()



  
		
           
