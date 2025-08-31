import os
import subprocess

# 指定物种树文件名
species_tree = "rax.tree"

# 让用户输入外类群
outgroup = "Tupaia_chinensis_GCF_000334495.pep"

# 遍历当前文件夹中的所有文件
for filename in os.listdir("."):
    if filename.endswith(".fasta"):  # 只处理FASTA文件
        # 从文件名中提取编号，例如 file1.fasta -> 1
        num = filename.split(".")[0][4:]  # 假设文件名格式为fileX.fasta
        # 构造RAxML命令，使用10个线程，并指定外类群
        raxml_cmd = f"raxmlHPC-PTHREADS -T 10 -f e -t {species_tree} -s {filename} -m GTRGAMMA -o {outgroup} -n output{num}"
        # 执行RAxML命令
        subprocess.run(raxml_cmd, shell=True)
        # 重命名RAxML生成的树文件
        result_tree = f"RAxML_result.output{num}"
        os.rename(result_tree, f"tree{num}.newick")

print("处理完成，已生成物种树！")