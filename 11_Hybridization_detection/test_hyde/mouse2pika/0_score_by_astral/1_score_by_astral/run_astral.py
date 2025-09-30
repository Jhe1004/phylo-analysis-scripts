import os

now_dir = os.getcwd()
input_trees_dir = "input"  # 输入文件所在的子目录

def get_topology_file():
    # 获取当前文件夹中指定拓扑结构的树文件
    # 输入：无
    # 输出：单个拓扑树文件路径（假设只有一个符合条件的文件）
    file_temp = os.listdir(now_dir + os.sep + input_trees_dir)
    for each in file_temp:
        if ".topology" in each:  # 假设拓扑文件以 .topology 结尾，可根据实际情况调整
            return each
    return None  # 如果没有找到符合条件的文件，返回 None

def get_gene_trees_file():
    # 获取当前文件夹中的基因树文件
    # 输入：无
    # 输出：基因树文件路径（假设只有一个 .trees 文件）
    file_temp = os.listdir(now_dir + os.sep + input_trees_dir)
    for each in file_temp:
        if ".trees" in each:
            return each
    return None

def main():
    # 获取拓扑树文件和基因树文件
    topology_file = get_topology_file()
    gene_trees_file = get_gene_trees_file()
    
    if not topology_file:
        print("Error: No topology file found in the input directory.")
        return
    if not gene_trees_file:
        print("Error: No gene trees file (.trees) found in the input directory.")
        return
    
    # 构造 ASTRAL 命令
    input_gene_trees = os.path.join(input_trees_dir, gene_trees_file)  # 基因树文件路径
    input_topology = os.path.join(input_trees_dir, topology_file)      # 拓扑树文件路径
    output_file = topology_file.replace(".topology", "_with_branches.tree")  # 输出文件名
    
    # 使用 -q 选项为已有拓扑计算枝长
    astral_cmd = (
        f"java -jar astral.jar -i {input_gene_trees} -q {input_topology} "
        f"-o {output_file}"
    )
    print(f"Running command: {astral_cmd}")
    os.system(astral_cmd)

if __name__ == "__main__":
    main()


  
		
           
