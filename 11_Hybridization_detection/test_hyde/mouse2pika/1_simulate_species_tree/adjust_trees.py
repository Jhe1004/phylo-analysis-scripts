import os
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from io import StringIO

# 读取被移动的分支和目标分支的叶节点列表
def read_clades(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        if len(lines) < 2:
            print("错误：clades.txt 文件中应包含两行，分别为被移动的分支和目标分支的叶节点列表。")
            return None, None
        moved_clade_leaves = lines[0].strip().split(',')
        target_clade_leaves = lines[1].strip().split(',')
    return moved_clade_leaves, target_clade_leaves

# 打印树中的所有克莱德（分支）
def print_all_clades(tree):
    print("树中的克莱德（分支）列表：")
    clade_counter = 1  # 克莱德编号
    printed_clades = set()  # 用于记录已打印的克莱德，防止重复
    for clade in tree.find_clades(order='preorder'):
        if clade.is_terminal():
            continue
        # 获取该克莱德下的所有叶节点名称
        leaf_names = [leaf.name for leaf in clade.get_terminals()]
        # 将叶节点名称排序，便于比较
        leaf_names_sorted = tuple(sorted(leaf_names))
        # 如果该克莱德的叶节点组合尚未打印过，则打印
        if leaf_names_sorted not in printed_clades:
            printed_clades.add(leaf_names_sorted)
            print(f"克莱德 {clade_counter}: 包含的叶节点: {', '.join(leaf_names)}")
            clade_counter += 1
    print("\n")

# 调整树的函数
def adjust_tree(tree, moved_clade_leaves, target_clade_leaves):
    # 找到要移动的分支
    clade_to_move = tree.common_ancestor(moved_clade_leaves)
    if clade_to_move is None:
        print("错误：无法找到被移动的分支。")
        return tree

    # 移除要移动的分支
    parent_of_moved = clade_to_move.parent
    if parent_of_moved is None:
        print("错误：被移动的分支是根节点，无法移动。")
        return tree
    # 保存被移动分支的分支长度
    moved_branch_length = clade_to_move.branch_length if clade_to_move.branch_length is not None else 0.0
    parent_of_moved.clades.remove(clade_to_move)
    print("已从原位置移除被移动的分支。")

    # 检查并删除冗余节点
    while True:
        if len(parent_of_moved.clades) == 0:
            # 如果父节点没有子节点，删除父节点
            grandparent = parent_of_moved.parent
            if grandparent is None:
                # 已经到达根节点
                print("已达到根节点，停止删除冗余节点。")
                break
            # 保存父节点的分支长度
            parent_branch_length = parent_of_moved.branch_length if parent_of_moved.branch_length is not None else 0.0
            grandparent.clades.remove(parent_of_moved)
            print("删除了一个没有子节点的冗余节点。")
            parent_of_moved = grandparent
        elif len(parent_of_moved.clades) == 1:
            # 如果父节点只有一个子节点，将唯一子节点提升
            only_child = parent_of_moved.clades[0]
            child_branch_length = only_child.branch_length if only_child.branch_length is not None else 0.0
            parent_branch_length = parent_of_moved.branch_length if parent_of_moved.branch_length is not None else 0.0
            # 合并分支长度
            only_child.branch_length = child_branch_length + parent_branch_length
            # 将唯一子节点替换父节点
            grandparent = parent_of_moved.parent
            if grandparent is None:
                # 父节点是根节点，更新树的根节点
                tree.root = only_child
                only_child.parent = None
                print("树的根节点已更新。")
                break
            else:
                grandparent.clades.remove(parent_of_moved)
                grandparent.clades.append(only_child)
                only_child.parent = grandparent
                print("提升了一个只有一个子节点的冗余节点。")
            parent_of_moved = grandparent
        else:
            # 父节点有多个子节点，不需要删除
            break

    # 找到目标分支的位置
    target_clade = tree.common_ancestor(target_clade_leaves)
    if target_clade is None:
        print("错误：无法找到目标分支。")
        return tree

    parent_of_target = target_clade.parent
    if parent_of_target is None:
        print("错误：目标分支是根节点，无法在其上进行操作。")
        return tree

    # 保存目标分支的分支长度
    target_branch_length = target_clade.branch_length if target_clade.branch_length is not None else 0.0
    parent_of_target.clades.remove(target_clade)
    print("已从父节点移除目标分支。")

    # 创建一个新的节点，将目标节点和要移动的分支作为其子节点
    new_clade = Clade()
    # 设置新节点的分支长度为目标节点原来的分支长度
    new_clade.branch_length = target_branch_length
    # 将目标节点和要移动的分支添加为新节点的子节点
    new_clade.clades.append(target_clade)
    new_clade.clades.append(clade_to_move)
    # 设置新节点为目标节点和被移动节点的父节点
    target_clade.parent = new_clade
    clade_to_move.parent = new_clade
    print("已创建新的节点，并添加目标分支和被移动的分支为子节点。")

    # 将新节点添加回父节点
    parent_of_target.clades.append(new_clade)
    new_clade.parent = parent_of_target
    print("已将新的节点添加回父节点。")

    return tree

# 主函数，批量处理 .tree 文件
def process_trees(clade_file):
    # 获取当前目录下的所有 .tree 文件
    tree_files = [f for f in os.listdir('.') if f.endswith('.tree')]
    if not tree_files:
        print("当前目录下没有找到 .tree 文件。")
        return
    for tree_file in tree_files:
        print(f"\n正在处理文件：{tree_file}")
        with open(tree_file, 'r') as f:
            newick_str = f.read()
        tree = Phylo.read(StringIO(newick_str), 'newick')
        # 为每个节点设置 parent 属性，方便操作
        for clade in tree.find_clades():
            for child in clade.clades:
                child.parent = clade
        # 打印树中的所有克莱德
        print_all_clades(tree)
        # 提示用户输入按 Enter 键继续
        input("请查看以上克莱德列表，更新 clades.txt 文件后按 Enter 键继续...")
        # 读取被移动的分支和目标分支的叶节点列表
        moved_clade_leaves, target_clade_leaves = read_clades(clade_file)
        if not moved_clade_leaves or not target_clade_leaves:
            print("错误：clades.txt 文件内容不完整。")
            continue
        # 调整树
        adjusted_tree = adjust_tree(tree, moved_clade_leaves, target_clade_leaves)
        # 将调整后的树保存为新的文件
        new_tree_file = tree_file.replace('.tree', '.new.tree')
        Phylo.write(adjusted_tree, new_tree_file, 'newick')
        print(f"已保存调整后的树到文件：{new_tree_file}")

if __name__ == "__main__":
    # 指定包含分支信息的文本文件
    clade_file = 'clades.txt'
    process_trees(clade_file)