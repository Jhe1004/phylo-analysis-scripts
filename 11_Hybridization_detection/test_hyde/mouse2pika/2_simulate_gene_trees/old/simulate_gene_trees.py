import subprocess
import os
import io
from Bio import Phylo

def convert_to_coalescent_units(tree, Ne, generation_time):
    """
    将物种树的分支长度从时间单位转换为Coalescent单位。
    """
    for clade in tree.find_clades():
        if clade.branch_length:  # 如果存在分支长度
            clade.branch_length = clade.branch_length / (2 * Ne * generation_time)
    return tree

def simulate_gene_trees(species_tree_file, ngene_trees, output_file):
    """
    使用coalescent模型调用R脚本模拟基因树。
    """
    with open(species_tree_file, 'r') as f:
        species_tree_str = f.read().strip()

    sp_tree = Phylo.read(io.StringIO(species_tree_str), "newick")
    sp_names = [leaf.name for leaf in sp_tree.get_terminals()]
    nspecies = len(sp_names)

    r_script = f"""
library(phybase)

args<-commandArgs(T)
con <- file(args[1], "r")
line=readLines(con)
close(con)
mptree1 = line

spname <- species.name(mptree1)
nodematrix <- read.tree.nodes(str=mptree1, name=spname)$nodes
nodematrix[,5]<-2

a3 <- as.numeric(args[2])
genetrees=1:a3
for(i in 1:a3)
genetrees[i]<-sim.coaltree.sp(rootnode=nrow(nodematrix), nodematrix=nodematrix, nspecies=length(spname), seq=rep(1,length(spname)), name=spname)$gt
write(genetrees, file="{output_file}")
"""

    with open("simulate_gene_trees.R", "w") as r_file:
        r_file.write(r_script)

    subprocess.run(["Rscript", "simulate_gene_trees.R", species_tree_file, str(ngene_trees)])
    os.remove("simulate_gene_trees.R")

def reroot_gene_trees(gene_trees_file, outgroup_species, output_file):
    """
    使用指定的外类群物种重新置根模拟的基因树。
    """
    with open(gene_trees_file, "r") as read_file:
        with open(output_file, "w") as write_file:
            for each_line in read_file:
                tree = Phylo.read(io.StringIO(each_line.strip()), "newick")
                try:
                    tree.root_with_outgroup(outgroup_species)
                    handle = io.StringIO()
                    Phylo.write(tree, handle, "newick")
                    tree_str = handle.getvalue().strip()
                    write_file.write(tree_str + "\n")
                except Exception as e:
                    print(f"无法在以下树中设置外类群：\n{each_line}\n错误信息：{e}")

def main():
    # 直接定义参数
    species_tree = "old.tree"  # 请改为你的物种树文件路径
    Ne = 10000  # 有效种群大小
    generation_time = 0.000002  # 世代长度（百万年）
    ngene_trees = 900  # 要模拟的基因树数量
    outgroup_species = "Tupaia_chinensis_GCF_000334495.pep"  # 请改为你的外类群物种名称

    # 从输入的物种树文件名生成输出文件前缀
    base_filename = os.path.splitext(os.path.basename(species_tree))[0]

    # 第一步：将分支长度转换为Coalescent单位
    tree = Phylo.read(species_tree, "newick")
    converted_tree = convert_to_coalescent_units(tree, Ne, generation_time)
    converted_tree_file = f"{base_filename}_coalescent_units_tree.newick"
    handle = io.StringIO()
    Phylo.write(converted_tree, handle, "newick")
    converted_tree_str = handle.getvalue().strip()
    with open(converted_tree_file, 'w') as f:
        f.write(converted_tree_str)
    print(f"转换后的物种树已保存到 {converted_tree_file}")

    # 第二步：使用R脚本模拟基因树
    simulate_gene_trees(converted_tree_file, ngene_trees, f"{base_filename}_simulate.trees")
    print(f"模拟的基因树已保存到 {base_filename}_simulate.trees")

    # 第三步：重新置根基因树
    reroot_gene_trees(f"{base_filename}_simulate.trees", outgroup_species, f"{base_filename}_rerooted_gene_trees.tree")
    print(f"重新置根的基因树已保存到 {base_filename}_rerooted_gene_trees.tree")

if __name__ == "__main__":
    main()