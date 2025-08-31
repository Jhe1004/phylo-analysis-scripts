from ete3 import Tree

# 定义参数（请根据你的实际情况修改路径）
INPUT_TREE_FILE = "pika_with_branches.tree"  # 输入树文件路径
OUTPUT_TREE_FILE = "pika_with_branches_u.tree"  # 输出树文件路径

def make_ultrametric(tree):
    """
    将树转换为超度量树，使所有叶子到根的距离相等。
    
    参数:
        tree (ete3.Tree): 输入的树对象
    
    返回:
        ete3.Tree: 转换后的超度量树
    """
    # 计算每个叶子到根的距离
    leaf_distances = {}
    for leaf in tree:
        distance = leaf.get_distance(tree)  # 获取到根的距离
        leaf_distances[leaf.name] = distance
    
    # 找到最大距离，作为目标距离
    max_distance = max(leaf_distances.values())
    
    # 调整每个叶子的枝长
    for leaf in tree:
        current_distance = leaf_distances[leaf.name]
        if current_distance < max_distance:
            # 计算需要延长的枝长
            extension = max_distance - current_distance
            # 延长叶子的枝长
            leaf.dist += extension
    
    return tree

def main():
    # 读取输入树
    tree = Tree(INPUT_TREE_FILE)
    
    # 转换为超度量树
    ultrametric_tree = make_ultrametric(tree)
    
    # 保存新的树到文件
    ultrametric_tree.write(outfile=OUTPUT_TREE_FILE, format=1)  # format=1 保留枝长
    
    print(f"超度量树已保存到 {OUTPUT_TREE_FILE}")

if __name__ == "__main__":
    main()