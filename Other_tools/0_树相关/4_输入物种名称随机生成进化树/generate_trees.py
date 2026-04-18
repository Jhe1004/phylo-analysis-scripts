#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
from pathlib import Path

from ete3 import Tree


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/0_树相关/4_输入物种名称随机生成进化树/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/0_树相关/4_输入物种名称随机生成进化树/output"
INPUT_TAXA_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/0_树相关/4_输入物种名称随机生成进化树/input/taxa.txt"
OUTGROUP_TAXON_NAME = "Clematis_songorica_Altay"
NUM_TREES_TO_GENERATE = 10


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = Path(OUTPUT_DIRECTORY) / "random_trees"


def main() -> None:
    taxa_file = INPUT_DIR / INPUT_TAXA_FILE
    if not taxa_file.exists():
        raise FileNotFoundError(f"未找到物种列表文件: {taxa_file}")
    taxon_names = [line.strip() for line in taxa_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    if OUTGROUP_TAXON_NAME not in taxon_names:
        raise ValueError(f"指定的外类群不在物种列表中: {OUTGROUP_TAXON_NAME}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for i in range(NUM_TREES_TO_GENERATE):
        t = Tree()
        t.populate(len(taxon_names), names_library=random.sample(taxon_names, len(taxon_names)), random_branches=True)
        t.resolve_polytomy()
        outgroup_leaf = t.get_leaves_by_name(OUTGROUP_TAXON_NAME)[0]
        t.set_outgroup(outgroup_leaf)
        t.write(outfile=str(OUTPUT_DIR / f"random_tree_{i+1}.nwk"), format=9)
    print("随机树生成完成。")


if __name__ == "__main__":
    main()
