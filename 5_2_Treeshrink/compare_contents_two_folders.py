'''
该脚本可以比较Treeshrink相关中两个文件夹中包含的树文件和phy文件是否严格的一一对应
'''


import os
now_dir = os.getcwd()

'''
得到目录combine_tree_with_name_file中的best_tree名称，并将其后面的数字识别号存入到变量best_tree_file_name中
'''

def get_tree_list():            #得到当前目录   
  file_temp = os.listdir(path=now_dir + "/combine_tree_with_name")    #得到当前目录中的所有文件名
  best_tree_file_name = []         #储存基因文件名称
  for each in file_temp:        #将后缀为.fasta的文件加入列表file_name中
      if "RAxML_bestTree" in each:
          best_tree_file_name.append(each[20:-12])
  return best_tree_file_name


'''
得到目录filtrate_alienment_use_treeshrink_output中的phy名称，并将其后面的数字识别号存入到变量phy_file_name中
'''

def get_alignment_list(): 
  file_temp = os.listdir(path=now_dir + "/filtrate_alignment_use_treeshrink_output")    #得到当前目录中的所有文件名
  phy_file_name = []         #储存基因文件名称
  for each in file_temp:        #将后缀为.fasta的文件加入列表file_name中
      if "cds_maffted.fasta" in each:
          phy_file_name.append(each[5:-18])
  return phy_file_name  




best_tree_file_list = get_tree_list()
alignment_file_list = get_alignment_list()

'''
得到两个列表，并且互相比较，将比较出的冲突的地方输出
'''

print("best_tree_file_list 数目 " + str(len(best_tree_file_list)))
print("alignment_file_list 数目 " + str(len(alignment_file_list)))


for each_best_tree_file in best_tree_file_list:
  if each_best_tree_file not in alignment_file_list:
    print("best_tree_file " + each_best_tree_file + " not in phy file list")
  
for each_phy_file in alignment_file_list:
  if each_phy_file not in best_tree_file_list:
    print("phy_file " + each_phy_file + " not in each_best_tree_file")
