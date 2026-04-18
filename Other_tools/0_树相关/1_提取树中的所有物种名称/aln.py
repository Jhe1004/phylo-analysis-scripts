#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from ete3 import Tree


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/0_树相关/1_提取树中的所有物种名称/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/0_树相关/1_提取树中的所有物种名称/output"
TREE_EXTENSIONS = [".tree", ".nwk", ".newick"]


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tree_files = []
    for ext in TREE_EXTENSIONS:
        tree_files.extend(sorted(INPUT_DIR.glob(f"*{ext}")))
    if not tree_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到树文件。")
    for input_tree in tree_files:
        t = Tree(str(input_tree))
        output_txt = OUTPUT_DIR / f"{input_tree.name}.txt"
        output_txt.write_text("\n".join(t.get_leaf_names()) + "\n", encoding="utf-8")
    print("树中物种名称提取完成。")


if __name__ == "__main__":
    main()
