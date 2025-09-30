'''
计算物种树各节点被多少基因树支持和不支持
基因树均存放在gene_tree文件夹中，每棵树一个文件
物种树存放在species_tree文件夹中，树的名称需要在脚本中修改
'''
import os

gene_tree_dir = "gene_tree"
species_tree = "species_tree/cds_min_seq_10.tree.newick"

os.system("java -jar target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 0 -d " +  gene_tree_dir + " -m " + species_tree + " -o out")