from pathlib import Path

from ete3 import Tree


# ==================== 配置区（直接修改这里）====================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/0_将想要检测的树提取出来/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/0_将想要检测的树提取出来/output"

INPUT_TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/0_将想要检测的树提取出来/input/result.tree"
SUBTREE_SPECIES_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/0_将想要检测的树提取出来/input/list.txt"
OUTGROUP_SPECIES = "Gossypium_hirsutum_SRR8156069.fasta.transdecoder.pep"
TREE_FORMAT = 0
# ============================================================


OUTPUT_TREE_FILE = "out.trees"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def read_species_list(species_file: Path) -> list[str]:
    return [line.strip() for line in species_file.read_text(encoding="utf-8").splitlines() if line.strip()]


def main() -> None:
    input_tree_path = INPUT_DIR / INPUT_TREE_FILE
    species_file_path = INPUT_DIR / SUBTREE_SPECIES_FILE
    output_tree_path = OUTPUT_DIR / OUTPUT_TREE_FILE

    if not input_tree_path.exists():
        raise FileNotFoundError(f"未找到输入树文件: {input_tree_path}")
    if not species_file_path.exists():
        raise FileNotFoundError(f"未找到待保留物种列表文件: {species_file_path}")

    subtree_species_list = read_species_list(species_file_path)
    if OUTGROUP_SPECIES not in subtree_species_list:
        raise ValueError("OUTGROUP_SPECIES 必须包含在待提取物种列表中。")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    kept_trees = []
    skipped_count = 0
    for raw_tree in input_tree_path.read_text(encoding="utf-8").splitlines():
        raw_tree = raw_tree.strip()
        if not raw_tree:
            continue
        try:
            tree = Tree(raw_tree, format=TREE_FORMAT)
            tree.prune(subtree_species_list, preserve_branch_length=True)
            tree.set_outgroup(OUTGROUP_SPECIES)
            kept_trees.append(tree.write(format=0).strip())
        except Exception:
            skipped_count += 1

    output_tree_path.write_text(
        "\n".join(kept_trees) + ("\n" if kept_trees else ""),
        encoding="utf-8",
    )

    print(f"提取完成，共保留 {len(kept_trees)} 棵子树，跳过 {skipped_count} 棵无法处理的树。")
    print(f"输出文件: {output_tree_path}")


if __name__ == "__main__":
    main()
