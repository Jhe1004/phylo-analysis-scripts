import os
import pyvolve
from io import StringIO
from Bio import Phylo

# 自动读取当前文件夹中的所有 .trees 文件
tree_files = [f for f in os.listdir('.') if f.endswith('.tree')]

# 定义GTR模型参数
gtr_params = {
    "kappa": 2.0,  # 转换/颠换比
    "state_freqs": [0.25, 0.25, 0.25, 0.25],  # A, C, G, T 的基频
    "rate_matrix": [
        [ "-", 1, 1, 1 ],
        [ 1, "-", 1, 1 ],
        [ 1, 1, "-", 1 ],
        [ 1, 1, 1, "-" ]
    ]
}

# 创建GTR模型
model = pyvolve.Model("nucleotide", params=gtr_params)

# 定义枝长缩放因子
scale_factor = 1

# 处理每个 .trees 文件
for tree_file in tree_files:
    # 读取基因树文件
    with open(tree_file, 'r') as f:
        tree_strings = f.read().strip().split(';')  # 按分号分割每棵树

    # 移除可能的空字符串
    tree_strings = [t.strip() for t in tree_strings if t.strip()]

    # 遍历该文件中的每棵基因树并进行序列模拟
    for idx, tree_str in enumerate(tree_strings):
        # 添加必要的分号
        if not tree_str.endswith(';'):
            tree_str += ';'

        # 使用 Bio.Phylo 解析树
        tree = Phylo.read(StringIO(tree_str), "newick")

        # 缩放枝长
        for clade in tree.find_clades():
            if clade.branch_length is not None:
                clade.branch_length *= scale_factor

        # 将缩放后的树转换为 Newick 字符串
        scaled_tree_str = tree.format("newick")

        # 使用 pyvolve 读取缩放后的树
        scaled_tree = pyvolve.read_tree(tree=scaled_tree_str)

        # 设置序列长度为 1000 个碱基
        seq_length = 1000000

        # 定义模拟参数
        partition = pyvolve.Partition(models=model, size=seq_length)

        # 进行模拟
        evolver = pyvolve.Evolver(partitions=partition, tree=scaled_tree)

        # 生成输出文件名
        output_filename = f'{os.path.splitext(tree_file)[0]}_simulated_sequences_{idx+1}.fasta'

        # 执行模拟并将结果写入文件
        evolver(seqfile=output_filename, seqformat='fasta')

        print(f'已为文件 {tree_file} 中的第 {idx+1} 棵树模拟序列，并保存到 {output_filename}')