from pathlib import Path


# ==================== 配置区（直接修改这里）====================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/1_Get_phylonet_input/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/1_Get_phylonet_input/output"

INPUT_TREES_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/11_Hybridization_detection/Phylonet/1_Get_phylonet_input/input/out.trees"
COMPUTING_METHOD = "ML"  # 可选：MP、ML、MPL
MAX_ALLOW_GENE_FLOW_NUM = 5
PHYLONET_OPTIONS = "-x 10 -pl 20 -di"
OUTPUT_PREFIX = "out"
# ============================================================


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY
VALID_METHODS = {"MP", "ML", "MPL"}


def build_nexus_text(tree_lines: list[str], gene_flow_num: int) -> str:
    tree_block = []
    tree_names = []
    for index, tree_line in enumerate(tree_lines, start=1):
        tree_name = f"geneTree{index}"
        tree_block.append(f"Tree {tree_name} = {tree_line}")
        tree_names.append(tree_name)

    infer_command = (
        f"InferNetwork_{COMPUTING_METHOD} "
        f"({','.join(tree_names)}) {gene_flow_num} {PHYLONET_OPTIONS};"
    )
    return (
        "#NEXUS\n\n"
        "BEGIN TREES;\n\n"
        + "\n".join(tree_block)
        + "\nEND;\n\n"
        "BEGIN PHYLONET;\n\n"
        + infer_command
        + "\n\nEND;\n"
    )


def main() -> None:
    if COMPUTING_METHOD not in VALID_METHODS:
        raise ValueError(f"COMPUTING_METHOD 只能是 {sorted(VALID_METHODS)} 之一。")

    input_trees_path = INPUT_DIR / INPUT_TREES_FILE
    if not input_trees_path.exists():
        raise FileNotFoundError(f"未找到输入树文件: {input_trees_path}")

    tree_lines = [line.strip() for line in input_trees_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not tree_lines:
        raise ValueError(f"输入树文件为空: {input_trees_path}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for gene_flow_num in range(MAX_ALLOW_GENE_FLOW_NUM + 1):
        output_path = OUTPUT_DIR / f"{OUTPUT_PREFIX}_max_{gene_flow_num}_conf.nex"
        output_path.write_text(build_nexus_text(tree_lines, gene_flow_num), encoding="utf-8")
        print(f"已生成: {output_path}")


if __name__ == "__main__":
    main()
