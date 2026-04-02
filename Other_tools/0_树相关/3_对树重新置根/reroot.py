#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from ete3 import Tree


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
TREE_SUFFIX = ".trees"
OUTGROUP_LIST = ["Clematis_gratopsis.fasta.transdecoder.pep"]


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tree_files = sorted(INPUT_DIR.glob(f"*{TREE_SUFFIX}"))
    if not tree_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 *{TREE_SUFFIX} 文件。")
    for input_file in tree_files:
        output_path = OUTPUT_DIR / f"{input_file.name}_reroot.tree"
        count = 0
        with input_file.open("r", encoding="utf-8") as read_file, output_path.open("w", encoding="utf-8") as write_file:
            for line in read_file:
                line = line.strip()
                if not line:
                    continue
                t = Tree(line)
                outgroup_nodes = []
                for name in OUTGROUP_LIST:
                    outgroup_nodes.extend(t.get_leaves_by_name(name))
                if not outgroup_nodes:
                    continue
                ancestor = t.get_common_ancestor(outgroup_nodes)
                t.set_outgroup(ancestor)
                write_file.write(t.write() + "\n")
                count += 1
        print(f"{input_file.name}: 成功置根 {count} 棵树。")


if __name__ == "__main__":
    main()
