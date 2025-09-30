'''
批量的生成treePL软件所需要的配置文件 (纯ete3版本)
输入：指定后缀的树文件（如.new, .tree），文件中含有一棵树，所有树需要含有完全相同的物种名称。
需要一个已经用MCMCTree等软件dating好的物种树作为校准来源。
输出：为每个输入的树文件生成一个treePL的配置文件。
'''
import os
from ete3 import Tree

# --- 用户配置 ---
# 请将此处的 "species.dated.tree" 替换为你的校准物种树文件名
# 这个文件应该是已经完成分子钟定年的、超度量（ultrametric）的树
calibrated_tree_file = "dna.tree"
# 请将此处的 ".new" 替换为你的基因树文件的后缀
tree_file_extension = ".new"
# -----------------

def get_file_list(extension):
    """获取当前目录中所有具有指定后缀的文件列表。"""
    file_list = [f for f in os.listdir() if f.endswith(extension)]
    print(f"找到了 {len(file_list)} 个 '{extension}' 后缀的树文件。")
    return file_list

def get_calibrat_age(calibrated_tree_obj, monophy):
    """
    从校准物种树(ete3 Tree object)中获取指定单系群的年龄。
    对于一个超度量树，一个节点的年龄等于该节点到其任何一个后代叶子节点的距离。
    """
    if not monophy:
        print("警告: 尝试为一个空的类群列表获取年龄，返回0。")
        return 0.0
    try:
        # 使用ete3的get_common_ancestor方法找到MRCA（最近公共祖先）
        mrca = calibrated_tree_obj.get_common_ancestor(monophy)
        
        # 获取该MRCA节点下的任意一个叶子节点
        # 在超度量树中，从此节点到它所有后代叶子的距离都相等
        any_leaf = mrca.get_leaves()[0]
        
        # 计算MRCA节点到这个叶子的距离，即为该节点的年龄
        age = mrca.get_distance(any_leaf)
        return age
    except Exception as e:
        # 如果物种名称在校准树中找不到，会触发异常
        print(f"警告: 无法为类群 {monophy} 计算年龄。可能原因：物种名不匹配。错误: {e}")
        return 0.0

def generate_calibra_list(each_tree, calibrated_tree_obj):
    """为单个基因树生成校准点列表。"""
    calibra_list = []
    try:
        t = Tree(each_tree, format=1) # format=1 表示灵活的Newick格式
    except Exception as e:
        print(f"错误：无法读取树文件 {each_tree}。请检查文件格式。错误: {e}")
        return []

    def check_monophyly(gene_tree, calibrated_tree, node):
        """检查基因树中的一个节点在物种树中是否也是单系的。"""
        ingroup_leaves = node.get_leaf_names()
        
        if len(ingroup_leaves) <= 1:
            return False
            
        outgroup_leaves = [leaf for leaf in gene_tree.get_leaf_names() if leaf not in ingroup_leaves]
        
        try:
            species_mrca = calibrated_tree.get_common_ancestor(ingroup_leaves)
        except Exception:
            return False

        for leaf in species_mrca.get_leaf_names():
            if leaf in outgroup_leaves:
                return False
        return True

    # 1. 校准根节点
    all_leaves = t.get_leaf_names()
    root_age = get_calibrat_age(calibrated_tree_obj, all_leaves)
    if root_age > 0:
        calibra_list.append(f"mrca = root_node {' '.join(all_leaves)}\n")
        calibra_list.append(f"fixage = root_node {root_age}\n")

    # 2. 校准其他内部节点
    node_counter = 1
    for each_node in t.traverse():
        if not each_node.is_leaf():
            if check_monophyly(t, calibrated_tree_obj, each_node):
                node_leaves = each_node.get_leaf_names()
                node_age = get_calibrat_age(calibrated_tree_obj, node_leaves)
                
                if node_age > 0:
                    node_name = f"node_{node_counter}"
                    calibra_list.append(f"mrca = {node_name} {' '.join(node_leaves)}\n")
                    calibra_list.append(f"min = {node_name} {node_age * 0.95}\n")
                    calibra_list.append(f"max = {node_name} {node_age * 1.05}\n")
                    node_counter += 1
    return calibra_list

def generate_treePL_conf(each_tree, calibrated_tree_obj):
    """为单个树文件生成完整的treePL配置文件。"""
    conf_filename = each_tree.replace(tree_file_extension, ".conf")
    
    fasta_filename = each_tree.replace(tree_file_extension, ".fasta")
    seq_len = 0
    try:
        with open(fasta_filename) as read_file:
            for each_line in read_file:
                if not each_line.startswith(">") and each_line.strip():
                    seq_len = len(each_line.strip())
                    break
    except FileNotFoundError:
        print(f"警告: 找不到对应的FASTA文件 '{fasta_filename}'。'numsites' 将被设置为0。")
        print("请手动在生成的配置文件中修改 'numsites' 的值。")

    with open(conf_filename, "w") as write_file:
        write_file.write(f"// Configuration file for {each_tree}\n\n")
        
        write_file.write("[Input files containing the ML trees]\n")
        write_file.write(f"treefile = {each_tree}\n\n")
        
        write_file.write("[General commands]\n")
        write_file.write("nthreads = 12\n")
        write_file.write("smooth = 1\n")
        write_file.write("thorough\n")
        write_file.write("log_pen\n")
        if seq_len > 0:
            write_file.write(f"numsites = {seq_len}\n\n")
        else:
            write_file.write(f"// numsites = 0 // WARNING: Please set sequence length manually!\n\n")

        write_file.write("[Calibrations]\n")
        calib_lines = generate_calibra_list(each_tree, calibrated_tree_obj)
        if not calib_lines:
            write_file.write("// No valid calibration points found.\n")
        else:
            for each_line in calib_lines:
                write_file.write(each_line)
        write_file.write("\n")

        write_file.write("[Optimisation parameters]\n")
        write_file.write("opt = 3\n")
        write_file.write("moredetail\n")
        write_file.write("optad = 3\n")
        write_file.write("moredetailad\n")
        write_file.write("optcvad = 5\n")
        write_file.write("randomcv\n")
        write_file.write("cviter = 5\n")
        write_file.write("cvsimaniter = 100000\n")
        write_file.write("cvstart = 10000\n")
        write_file.write("cvstop = 0.00001\n")
        write_file.write("cvmultstep = 0.1\n")
        write_file.write(f"cvoutfile = {each_tree}.cvout\n")
        write_file.write(f"outfile = {each_tree}.newtree\n")
        
    print(f"已生成配置文件: {conf_filename}")

def main():
    """主函数，执行整个流程。"""
    tree_files = get_file_list(tree_file_extension)
    if not tree_files:
        print(f"在当前目录下未找到任何 '{tree_file_extension}' 文件，程序退出。")
        return
        
    try:
        print(f"正在加载校准物种树: {calibrated_tree_file}")
        # 在主函数中只加载一次校准树，提高效率
        calibrated_tree_obj = Tree(calibrated_tree_file, format=1)
        print("校准物种树加载成功。")
    except Exception as e:
        print(f"致命错误: 无法加载校准物种树 '{calibrated_tree_file}'。")
        print(f"请检查文件名是否正确，文件是否存在且为有效的Newick格式。错误详情: {e}")
        return
        
    for each_tree in tree_files:
        print(f"\n--- 正在处理: {each_tree} ---")
        generate_treePL_conf(each_tree, calibrated_tree_obj)

if __name__ == "__main__":
    main()
