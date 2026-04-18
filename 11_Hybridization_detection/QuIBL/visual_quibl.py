import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ete3 import Tree


# ==================== 配置区（直接修改这里）====================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/QuIBL/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/QuIBL/output"

TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/QuIBL/input/concat.newick"
QUIBL_RESULT_CSV = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/QuIBL/input/quibl_results.csv"
LOG_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/QuIBL/output/result_log.txt"
HEATMAP_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/QuIBL/output/heatmap_with_tree.png"
BIC_DIFFERENCE_THRESHOLD = 20.0
COLORMAP = "hot"
# ============================================================


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def parse_triplet(triplet_str: str, leaf_names: list[str]) -> list[str]:
    tokens = triplet_str.split("_")
    result = []
    start = 0
    while start < len(tokens):
        matched = False
        for end in range(start + 1, len(tokens) + 1):
            candidate = "_".join(tokens[start:end])
            if candidate in leaf_names:
                result.append(candidate)
                start = end
                matched = True
                break
        if not matched:
            raise ValueError(f"无法解析 triplet: {triplet_str}")

    if len(result) != 3:
        raise ValueError(f"triplet 解析结果不是 3 个物种: {triplet_str} -> {result}")
    return result


def is_outgroup_correct(tree: Tree, sp1: str, sp2: str, outgroup: str) -> bool:
    pruned_tree = tree.copy()
    pruned_tree.prune([sp1, sp2, outgroup], preserve_branch_length=True)
    expected_tree = Tree(f"(({sp1},{sp2}),{outgroup});")
    expected_tree.set_outgroup(expected_tree.search_nodes(name=outgroup)[0])
    return pruned_tree.robinson_foulds(expected_tree)[0] == 0


def main() -> None:
    tree_path = INPUT_DIR / TREE_FILE
    quibl_result_path = INPUT_DIR / QUIBL_RESULT_CSV
    if not tree_path.exists():
        raise FileNotFoundError(f"未找到物种树文件: {tree_path}")
    if not quibl_result_path.exists():
        raise FileNotFoundError(f"未找到 QuIBL 结果文件: {quibl_result_path}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    log_path = OUTPUT_DIR / LOG_FILE
    heatmap_path = OUTPUT_DIR / HEATMAP_FILE

    species_tree = Tree(tree_path.read_text(encoding="utf-8").strip(), format=1)
    leaf_names = species_tree.get_leaf_names()
    result_dict: dict[str, list[float]] = {}

    with quibl_result_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            if not row:
                continue
            triplet = row[0]
            outgroup = row[1]
            species_triplet = parse_triplet(triplet, leaf_names)
            pair = [species for species in species_triplet if species != outgroup]
            if len(pair) != 2:
                continue

            sp1, sp2 = pair
            if is_outgroup_correct(species_tree, sp1, sp2, outgroup):
                continue

            bic2_dist = float(row[8])
            bic1_dist = float(row[9])
            diff = bic1_dist - bic2_dist
            mixprop2 = float(row[5])
            value = mixprop2 if diff > BIC_DIFFERENCE_THRESHOLD else 0.0
            key = f"{sp1}:{sp2}"
            result_dict.setdefault(key, []).append(value)

    with log_path.open("w", encoding="utf-8") as log_handle:
        for key, values in sorted(result_dict.items()):
            average = sum(values) / len(values) if values else 0.0
            log_handle.write(f"{key}: {values}, Average: {average}\n")

    matrix_size = len(leaf_names)
    average_matrix = np.zeros((matrix_size, matrix_size))
    for row_index, species_1 in enumerate(leaf_names):
        for col_index, species_2 in enumerate(leaf_names):
            if row_index == col_index:
                continue
            key = (
                f"{species_1}:{species_2}"
                if f"{species_1}:{species_2}" in result_dict
                else f"{species_2}:{species_1}"
            )
            values = result_dict.get(key, [])
            average_matrix[row_index, col_index] = sum(values) / len(values) if values else 0.0

    figure, axis = plt.subplots(figsize=(12, 10))
    image = axis.imshow(average_matrix, cmap=COLORMAP, interpolation="nearest")
    figure.colorbar(image, ax=axis, label="Average Mixprop2")
    axis.set_xticks(range(matrix_size))
    axis.set_yticks(range(matrix_size))
    axis.set_xticklabels(leaf_names, rotation=90)
    axis.set_yticklabels(leaf_names)
    axis.set_title("Heatmap of Average Mixprop2")
    plt.tight_layout()
    plt.savefig(heatmap_path, dpi=300)
    plt.close(figure)

    print(f"可视化完成，热图文件: {heatmap_path}")
    print(f"日志文件: {log_path}")


if __name__ == "__main__":
    main()
