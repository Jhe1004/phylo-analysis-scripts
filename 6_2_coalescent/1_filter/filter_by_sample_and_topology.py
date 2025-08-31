'''
筛选出与参考基因树拓扑结构相同的基因树，所包含的物种名可以任意指定。
输入文件有alignment排序文件, 以及RAxML_bipartitions.*树文件，这两种文件需要放置在本文件夹中
'''
import os
from Bio import SeqIO
from ete3 import Tree
now_dir = os.getcwd()

reference_tree = "RAxML_bipartitions.result.newick"
outgoup = "Anemoclema_glauciifolium_19.fasta.transdecoder.pep"
samples =  "samples.txt"

ingroups = []
with open(samples, "r") as read_file:
  for each_line in read_file:
    ingroups.append(each_line.replace("\n", ""))


def get_file_list():
  #函数get_file_list：获取当前文件夹中指定文件
  #输入：无
  #输出：指定文件列表：file_list
  file_temp = os.listdir()
  file_list = []
  for each in file_temp:       
    if "_bipartitions.ortho" in each:
      file_list.append(each)
  return file_list

def get_subtree_in_reference():
  #根据samples中的名称, 将reference树中不属于这些名称的样品剔除
  tree = Tree(reference_tree)
  tree.prune(ingroups)
  return tree

def get_subtree_in_genetree(each_tree_file):
  #根据samples中的名称, 将基因树中不属于这些名称的样品剔除
  tree = Tree(each_tree_file)
  try:
    tree.prune(ingroups)
    #print(each_tree_file + " have all the samples")
    tree.set_outgroup(outgoup)
    return tree
  except:
    pass

def filter_sub_gene_tree():
  for each_tree_file in get_file_list():
    gene_subtree = get_subtree_in_genetree(each_tree_file)
    reference_subtree = get_subtree_in_reference()
    if gene_subtree:
      #print(gene_subtree)
      #print(reference_subtree)
      rf = gene_subtree.robinson_foulds(reference_subtree)
      if rf[0] < 3:
        print(each_tree_file)


filter_sub_gene_tree()


