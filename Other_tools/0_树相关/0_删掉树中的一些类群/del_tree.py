#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import shutil
from pathlib import Path

from ete3 import Tree


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
TREE_FILE_PREFIX = "dating."
DELETE_LIST_FILE = "tree.txt"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    delete_list_file = INPUT_DIR / DELETE_LIST_FILE
    if not delete_list_file.exists():
        raise FileNotFoundError(f"未找到删除名单文件: {delete_list_file}")
    del_leaves_set = {line.strip() for line in delete_list_file.read_text(encoding="utf-8").splitlines() if line.strip()}
    tree_files = sorted([p for p in INPUT_DIR.iterdir() if p.is_file() and p.name.startswith(TREE_FILE_PREFIX)])
    if not tree_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到以 {TREE_FILE_PREFIX} 开头的树文件。")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for tree_file in tree_files:
        output_file = OUTPUT_DIR / tree_file.name
        shutil.copy2(tree_file, output_file)
        t = Tree(str(tree_file), format=1)
        original_leaves = set(t.get_leaf_names())
        leaves_to_keep = [name for name in original_leaves if name not in del_leaves_set]
        if leaves_to_keep and len(leaves_to_keep) < len(original_leaves):
            t.prune(leaves_to_keep)
            t.write(format=1, outfile=str(output_file))
    print("删树类群完成。")


if __name__ == "__main__":
    main()
