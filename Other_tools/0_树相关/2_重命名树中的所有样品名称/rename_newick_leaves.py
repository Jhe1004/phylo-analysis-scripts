#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path

from ete3 import Tree


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/0_树相关/2_重命名树中的所有样品名称/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/0_树相关/2_重命名树中的所有样品名称/output"
TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/0_树相关/2_重命名树中的所有样品名称/input/input_tree.newick"
MAPPING_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/0_树相关/2_重命名树中的所有样品名称/input/mapping.csv"


OUTPUT_FILE = "renamed_tree.newick"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def create_replacement_map(mapping_file: Path) -> dict[str, str]:
    replacement_dict: dict[str, str] = {}
    with mapping_file.open("r", encoding="utf-8") as infile:
        reader = csv.reader(infile)
        for rows in reader:
            if len(rows) >= 2:
                replacement_dict[rows[0].strip()] = rows[1].strip()
    return replacement_dict


def main() -> None:
    tree_file = INPUT_DIR / TREE_FILE
    mapping_file = INPUT_DIR / MAPPING_FILE
    if not tree_file.exists():
        raise FileNotFoundError(f"未找到树文件: {tree_file}")
    if not mapping_file.exists():
        raise FileNotFoundError(f"未找到映射文件: {mapping_file}")
    name_map = create_replacement_map(mapping_file)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / OUTPUT_FILE

    replaced_count = 0
    total_leaves = 0
    with tree_file.open("r", encoding="utf-8") as infile, output_file.open("w", encoding="utf-8") as outfile:
        for newick_string in infile:
            newick_string = newick_string.strip()
            if not newick_string:
                continue
            tree = Tree(newick_string, format=1)
            for leaf in tree.iter_leaves():
                total_leaves += 1
                if leaf.name in name_map:
                    leaf.name = name_map[leaf.name]
                    replaced_count += 1
            outfile.write(tree.write(format=1) + "\n")
    print(f"完成。共检查 {total_leaves} 个叶子，替换 {replaced_count} 个名称。")


if __name__ == "__main__":
    main()
