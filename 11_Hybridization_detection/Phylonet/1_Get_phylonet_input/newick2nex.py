import os
from Bio import SeqIO
now_dir = os.getcwd()

input_trees = "out.trees"
max_allow_gene_flow_num = 5  #最大允许的基因流数目
computing_method = "ML" #推算系统发育网络所使用的方法，有MP、ML、MPL可选


def make_PhyloNet_nex_file():
  '''
  函数make_PhyloNet_nex_file: 根据.trees文件生成相应的PhyloNet输入文件
  输入：无
  输出：PhyloNet软件需要的Nex格式的输入文件
  '''
for each_max_num in range(max_allow_gene_flow_num + 1):
  with open(input_trees, "r") as read_file:
    n = 0
    first_block_list = []
    second_block_list = []
    for each_line in read_file:
      if len(each_line) > 1:
        n = n + 1
        first_block_list.append("Tree geneTree" + str(n) + " = " + each_line)
        second_block_list.append("geneTree" + str(n))
    with open(input_trees[:-6] + "_max_" + str(each_max_num) +  "_conf.nex", "a") as write_file:
      write_file.write("#NEXUS\n\nBEGIN TREES;\n\n")
      for each_first_block_tree in first_block_list:
        write_file.write(each_first_block_tree)
      write_file.write("END;\n\n\nBEGIN PHYLONET;\n\n"
                        "InferNetwork_" + computing_method + " (" + ",".join(second_block_list) + ") " + str(each_max_num) + " -x 10 -pl 20 -di;\n\nEND;")

make_PhyloNet_nex_file()




    