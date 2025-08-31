import os
from Bio import Phylo

# 获取当前文件夹中所有后缀为 .tree 的文件
tree_files = [f for f in os.listdir('.') if f.endswith('.tree')]

# 遍历每个 .tree 文件
for tree_file in tree_files:
    # 生成对应的输出文件名（后缀为 .trees）
    output_file = f"{os.path.splitext(tree_file)[0]}.trees"
    
    # 打开输出文件
    with open(output_file, 'w') as output_handle:
        # 读取输入文件中的所有树
        trees = Phylo.parse(tree_file, 'newick')
        
        # 遍历每棵树
        for tree in trees:
            # 遍历树中的每个分支，并将分支长度设置为 0.01
            for clade in tree.find_clades():
                if clade.branch_length is not None:  # 如果分支长度存在，则修改
                    clade.branch_length = 0.02
            
            # 将修改后的树写入输出文件
            Phylo.write(tree, output_handle, 'newick')
    
    print(f"已将 {tree_file} 中的分支长度修改为 0.0001，并保存至 {output_file}")

