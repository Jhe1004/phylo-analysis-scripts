#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from ete3 import Tree


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INGROUPS_FILE = "renamed.tree.txt"
TREE_NAME_PATTERN = "RAx"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    ingroups_file = INPUT_DIR / INGROUPS_FILE
    if not ingroups_file.exists():
        raise FileNotFoundError(f"未找到保留名单文件: {ingroups_file}")
    ingroups_list = [line.strip() for line in ingroups_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    tree_files = sorted([p for p in INPUT_DIR.iterdir() if p.is_file() and TREE_NAME_PATTERN in p.name])
    if not tree_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到包含 {TREE_NAME_PATTERN} 的树文件。")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for tree_file in tree_files:
        t = Tree(str(tree_file))
        keep_names = [name for name in t.get_leaf_names() if name in ingroups_list]
        if not keep_names:
            continue
        t.prune(keep_names)
        t.write(format=1, outfile=str(OUTPUT_DIR / f"{tree_file.name}_new.tree"))
    print("保留类群子树提取完成。")


if __name__ == "__main__":
    main()
