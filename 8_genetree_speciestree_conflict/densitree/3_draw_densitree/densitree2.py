'''
绘制云图（DensiTree），云图中的半透明细线为基因树集，
实线为物种树。此脚本解决了绘图时叶尖末端不齐的问题。
'''
import toytree
import toyplot
import toyplot.svg

# --- 用户配置 ---
# 包含多个基因树的文件 (例如treePL的输出合并文件或BEAST/MrBayes的.trees文件)
GENE_TREES_FILE = "genetree2.trees"
# 单个物种树文件 (例如MCMCTree或ASTRAL的输出)
SPECIES_TREE_FILE = "dna.tree"
# 输出的图片文件名
OUTPUT_SVG_FILE = "densi_tree_plot.svg"
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
print(f"正在加载物种树: {SPECIES_TREE_FILE}")
species_tree = toytree.tree(SPECIES_TREE_FILE)
print(f"正在加载基因树集: {GENE_TREES_FILE}")
gene_trees = toytree.mtree(GENE_TREES_FILE)

# 2. 确定叶节点的统一顺序
# 使用物种树的叶节点顺序作为标准，确保所有树的Y轴顺序一致
tip_order = species_tree.get_tip_labels()
print(f"已确定统一的叶节点顺序，共 {len(tip_order)} 个物种。")

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
    fixed_order=tip_order,
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
    fixed_order=tip_order,
    edge_type='c',      # 'c' for squared edges
    edge_style={
        "stroke": "black", # 物种树使用醒目的黑色
        "stroke-width": 1.0,
    },
    tip_labels_align=True, # 对齐叶节点标签
);

# 6. 保存图像
print(f"正在保存图像到: {OUTPUT_SVG_FILE}")
toyplot.svg.render(canvas, OUTPUT_SVG_FILE)
print("完成！")
