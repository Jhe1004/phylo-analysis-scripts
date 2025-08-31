from Bio import Phylo

def convert_to_coalescent_units(tree, Ne, generation_time):
    """
    将物种树的分支长度从时间单位转换为Coalescent单位。
    
    参数:
    tree: Biopython解析的物种树对象
    Ne: 有效种群大小
    generation_time: 世代长度
    
    返回:
    转换后的物种树
    """
    # 遍历树的每个分支并转换其长度
    for clade in tree.find_clades():
        if clade.branch_length:  # 如果存在分支长度
            # 使用公式: L_coal = T_time / (2 * Ne * generation_time)
            clade.branch_length = clade.branch_length / (2 * Ne * generation_time)
    return tree

def main():
    # 读取Newick格式的物种树文件
    input_tree_file = "species_tree.newick"
    output_tree_file = "coalescent_units_tree.newick"
    
    # 设定有效种群大小和世代长度
    Ne = 1000000  # 例如：有效种群大小为10,000
    generation_time = 0.000005  # 例如：每代是5年

    # 使用Biopython读取树
    tree = Phylo.read(input_tree_file, "newick")
    
    # 将分支长度转换为Coalescent单位
    converted_tree = convert_to_coalescent_units(tree, Ne, generation_time)
    
    # 输出转换后的树到新的Newick文件
    Phylo.write(converted_tree, output_tree_file, "newick")
    print(f"Converted tree saved to {output_tree_file}")

if __name__ == "__main__":
    main()