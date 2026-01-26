from ete3 import Tree, TreeStyle
import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec

tree_file = "concat.newick"  # 输入的Newick格式的物种树文件，用于构建和分析物种树
input_file = "testOut.csv"  # 输入的CSV文件，其中包含处理的三元组及相关信息
log_file = "result_log.txt"  # 结果日志文件，用于记录计算结果和平均值
heatmap_file = "heatmap_with_tree.png"  # 生成的热图文件

with open(tree_file, "r") as f:
    species_tree_str = f.read().strip()  # 去除首尾多余空白符

t = Tree(species_tree_str, format=1)
# 不再对 leaf_names sort，以树的实际叶顺序为准
leaf_names = t.get_leaf_names()

def parse_triplet(triplet_str, leaf_names):
    tokens = triplet_str.split("_")
    result = []
    start = 0
    while start < len(tokens):
        matched = False
        for end in range(start+1, len(tokens)+1):
            candidate = "_".join(tokens[start:end])
            if candidate in leaf_names:
                result.append(candidate)
                start = end
                matched = True
                break
        if not matched:
            raise ValueError(f"Cannot parse triplet: {triplet_str}")

    if len(result) != 3:
        raise ValueError(f"Not exactly three species parsed from {triplet_str}, got: {result}")
    return result

def is_outgroup_correct(tree, sp1, sp2, outg):
    pruned_tree = tree.copy()
    pruned_tree.prune([sp1, sp2, outg], preserve_branch_length=True)
    expected_tree = Tree(f"(({sp1},{sp2}),{outg});")
    expected_tree.set_outgroup(expected_tree.search_nodes(name=outg)[0])
    return pruned_tree.robinson_foulds(expected_tree)[0] == 0

result_dict = {}

with open(input_file, newline='') as f:
    reader = csv.reader(f)
    header = next(reader)

    for row in reader:
        triplet = row[0]
        outgroup = row[1]
        species_triplet = parse_triplet(triplet, leaf_names)
        pair = [sp for sp in species_triplet if sp != outgroup]
        if len(pair) != 2:
            continue
        sp1, sp2 = pair

        if is_outgroup_correct(t, sp1, sp2, outgroup):
            continue
        else:
            BIC2Dist = float(row[8])
            BIC1Dist = float(row[9])
            diff = BIC1Dist - BIC2Dist
            mixprop2 = float(row[5])
            val = mixprop2 if diff > 20 else 0
            key = f"{sp1}:{sp2}"
            if key not in result_dict:
                result_dict[key] = []
            result_dict[key].append(val)

# 将结果写入日志文件
with open(log_file, "w") as log:
    for key, values in result_dict.items():
        avg_value = sum(values) / len(values) if values else 0
        log.write(f"{key}: {values}, Average: {avg_value}\n")

# 根据物种树的顺序生成平均值矩阵
n = len(leaf_names)
avg_matrix = np.zeros((n, n))
for i, sp1 in enumerate(leaf_names):
    for j, sp2 in enumerate(leaf_names):
        if i == j:
            avg_matrix[i, j] = 0
        else:
            key = f"{sp1}:{sp2}" if f"{sp1}:{sp2}" in result_dict else f"{sp2}:{sp1}"
            values = result_dict.get(key, [])
            avg_matrix[i, j] = sum(values) / len(values) if values else 0

# 创建绘图布局
fig = plt.figure(figsize=(12, 10))
gs = gridspec.GridSpec(2, 2, width_ratios=[1,4], height_ratios=[1,4])


# 热图
ax_heatmap = fig.add_subplot(gs[1,1])
im = ax_heatmap.imshow(avg_matrix, cmap="hot", interpolation="nearest")
cb = fig.colorbar(im, ax=ax_heatmap, label="Average Mixprop2")

ax_heatmap.set_xticks(range(n))
ax_heatmap.set_yticks(range(n))
ax_heatmap.set_xticklabels(leaf_names, rotation=90)
ax_heatmap.set_yticklabels(leaf_names)

plt.suptitle("Heatmap of Average Mixprop2", fontsize=16)
plt.tight_layout()
plt.savefig(heatmap_file)
plt.show()