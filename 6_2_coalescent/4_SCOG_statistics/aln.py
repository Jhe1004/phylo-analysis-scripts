#!/usr/bin/env python3
import os
import re
from Bio import SeqIO
import dendropy

def get_gene_length(fasta_file):
    """计算FASTA文件中基因的长度"""
    total_length = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_length += len(record.seq)
    return total_length

def get_average_bootstrap(tree_file):
    """使用DendroPy计算系统发育树的平均bootstrap支持率"""
    tree = dendropy.Tree.get(path=tree_file, schema='newick', rooting='default-rooted', preserve_underscores=True)
    bootstrap_values = []

    for node in tree.nodes():
        if node.label and node.is_internal():
            try:
                bootstrap_value = float(node.label)
                bootstrap_values.append(bootstrap_value)
            except ValueError:
                pass  # 如果标签不是数字，忽略

    if bootstrap_values:
        average_bootstrap = sum(bootstrap_values) / len(bootstrap_values)
        return average_bootstrap
    else:
        return None

def main():
    fasta_pattern = re.compile(r'ortho\d+_cds_maffted\.fasta')
    tree_pattern = re.compile(r'RAxML_bipartitions\.ortho\d+_cds_maffted')

    results = []

    for filename in os.listdir('.'):
        if fasta_pattern.match(filename):
            gene_id = filename.split('_')[0]  # 提取基因编号，例如'ortho64'
            fasta_file = filename
            tree_file = f'RAxML_bipartitions.{gene_id}_cds_maffted'

            if os.path.exists(tree_file):
                gene_length = get_gene_length(fasta_file)
                average_bootstrap = get_average_bootstrap(tree_file)

                if average_bootstrap is not None:
                    results.append((gene_id, gene_length, average_bootstrap))
                else:
                    print(f"{gene_id}: 未找到bootstrap支持率信息。")
            else:
                print(f"{gene_id}: 对应的树文件{tree_file}不存在。")

    # 输出结果
    print("基因编号\t基因长度\t平均Bootstrap支持率")
    for gene_id, gene_length, average_bootstrap in results:
        print(f"{gene_id}\t{gene_length}\t{average_bootstrap:.2f}")

if __name__ == '__main__':
    main()