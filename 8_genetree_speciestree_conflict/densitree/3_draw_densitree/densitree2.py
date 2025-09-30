# -*- coding: utf-8 -*-

'''
绘制云图（DensiTree），云图中的半透明细线为基因树集，
实线为物种树。此脚本解决了绘图时叶尖末端不齐的问题，
并且允许用户通过外部文件指定叶节点的显示顺序。
'''
import toytree
import toyplot
import toyplot.svg
import sys # 导入sys模块以便在出错时退出

# --- 用户配置 ---
# 包含多个基因树的文件 (例如treePL的输出合并文件或BEAST/MrBayes的.trees文件)
GENE_TREES_FILE = "genetree.trees"
# 单个物种树文件 (例如MCMCTree或ASTRAL的输出)
SPECIES_TREE_FILE = "dated_rotated.tree"
# 用于指定叶节点顺序的文本文件 (每行一个叶节点名称)
TIP_ORDER_FILE = "dated_rotated.tree.txt"
# 输出的图片文件名
OUTPUT_SVG_FILE = "densi_tree_plot_correct_order.svg"
# -----------------

# --- 核心修正函数 ---
def scale_tree_to_height(tree, target_height=1.0):
    """
    手动将 toytree 对象的总高度缩放到一个目标值。
    这是一个健壮的、不受 toytree 版本更新影响的归一化方法。
    """
    # 获取树的根节点
    root = tree.treenode
    # 计算当前树的总高度
    current_height = root.height
    
    # 如果树没有高度，则直接返回，避免除以零
    if current_height == 0:
        return tree
        
    # 计算缩放因子
    scaling_factor = target_height / current_height
    
    # 遍历树中的每一个节点，将其分支长度乘以缩放因子
    for node in root.traverse():
        node.dist *= scaling_factor
        
    # toytree对象被直接修改，但返回它是一个好习惯
    return tree

# 1. 加载物种树和基因树
print("正在加载物种树: {}".format(SPECIES_TREE_FILE))
species_tree = toytree.tree(SPECIES_TREE_FILE)
print("正在加载基因树集: {}".format(GENE_TREES_FILE))
gene_trees = toytree.mtree(GENE_TREES_FILE)

# 2. 从外部文件加载叶节点的统一顺序
print("正在从文件加载指定的叶节点顺序: {}".format(TIP_ORDER_FILE))
try:
    with open(TIP_ORDER_FILE, 'r') as f:
        # 读取文件中的每一行，去除首尾空白字符，并过滤掉空行
        tip_order_from_file = [line.strip() for line in f if line.strip()]
    
    # 【关键修改】反转列表顺序
    # toytree绘图时，列表第一个元素在底部，最后一个在顶部。
    # 为了让文件第一行显示在图的顶部，需要将整个列表反转。
    tip_order = tip_order_from_file[::-1]
    
    print("已加载并反转指定的叶节点顺序，共 {} 个物种。".format(len(tip_order)))

except FileNotFoundError:
    print("错误: 找不到指定的顺序文件 '{}'。请检查文件名和路径。".format(TIP_ORDER_FILE))
    sys.exit(1) # 退出脚本

# 3. (关键步骤) 归一化所有树的深度
# 使用我们自己编写的函数来归一化树高
print("正在归一化所有树的深度...")
species_tree = scale_tree_to_height(species_tree, 1.0)
gene_trees.treelist = [scale_tree_to_height(t, 1.0) for t in gene_trees]
print("归一化完成。")

# 4. 创建画布和唯一的坐标轴
# 只需要一个坐标轴，所有图都画在上面
canvas = toyplot.Canvas(width=1000, height=800)
axes = canvas.cartesian(xlabel="Relative time (root to tip)")

# 5. 绘制云图和物种树
# a. 先绘制基因树云图
print("正在绘制基因树云图...")
gene_trees.draw_cloud_tree(
    axes=axes,
    fixed_order=tip_order, # 使用从文件加载并反转后的顺序
    edge_style={
        "stroke": toytree.colors[1],  # 使用柔和的颜色
        "stroke-opacity": 0.1,       # 设置透明度
        "stroke-width": 1.0,
    },
)

# b. 在云图之上绘制物种树
print("正在叠加绘制物种树...")
species_tree.draw(
    axes=axes,
    fixed_order=tip_order, # 同样使用从文件加载并反转后的顺序
    edge_type='c',         # 'c' for squared edges
    edge_style={
        "stroke": "black", # 物种树使用醒目的黑色
        "stroke-width": 1.0,
    },
    tip_labels_align=True, # 对齐叶节点标签
);

# 6. 保存图像
print("正在保存图像到: {}".format(OUTPUT_SVG_FILE))
toyplot.svg.render(canvas, OUTPUT_SVG_FILE)
print("完成！")