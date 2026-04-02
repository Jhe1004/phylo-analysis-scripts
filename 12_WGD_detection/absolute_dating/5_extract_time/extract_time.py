#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
from collections import defaultdict

import dendropy
from ete3 import Tree


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

INPUT_EXTENSION = ".newtree"
OUTPUT_FILE_SUFFIX = ".pep"
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def get_species_in_clade(node) -> list[str]:
    result = []
    for leaf in node.get_leaf_names():
        species_name = leaf.split("++")[0]
        if species_name not in result:
            result.append(species_name)
    return result


def get_node_age(tree_path: Path, monophyletic_group: list[str]) -> float:
    tree = dendropy.Tree.get_from_path(str(tree_path), "newick", preserve_underscores=True)
    tree.calc_node_ages()
    mrca = tree.mrca(taxon_labels=monophyletic_group)
    return float(mrca.age)


def extract_ages_from_tree(tree_path: Path) -> dict[str, list[float]]:
    ete_tree = Tree(str(tree_path))
    species_list = get_species_in_clade(ete_tree)
    result: dict[str, list[float]] = defaultdict(list)
    for species_name in species_list:
        for node in ete_tree.traverse():
            if node.is_leaf() or len(node.children) != 2:
                continue
            child0_species = get_species_in_clade(node.children[0])
            child1_species = get_species_in_clade(node.children[1])
            if species_name in child0_species and species_name in child1_species:
                result[species_name].append(get_node_age(tree_path, node.get_leaf_names()))
    return result


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    ages_dir = output_dir / "sample_ages"
    ensure_directory(ages_dir)

    tree_files = sorted(input_dir.glob(f"*{INPUT_EXTENSION}"))
    if not tree_files:
        raise FileNotFoundError(f"未在 {input_dir} 中找到 {INPUT_EXTENSION} 文件")

    merged: dict[str, list[float]] = defaultdict(list)
    for tree_file in tree_files:
        tree_ages = extract_ages_from_tree(tree_file)
        for species_name, ages in tree_ages.items():
            merged[species_name].extend(ages)

    summary_file = output_dir / "extract_time_summary.tsv"
    with open(summary_file, "w", encoding="utf-8") as summary_handle:
        summary_handle.write("sample\tcount\n")
        for species_name, ages in sorted(merged.items()):
            output_file = ages_dir / f"{species_name}{OUTPUT_FILE_SUFFIX}"
            with open(output_file, "w", encoding="utf-8") as handle:
                for age in ages:
                    handle.write(f"{age}\n")
            summary_handle.write(f"{species_name}\t{len(ages)}\n")

    print(f"输出样本数: {len(merged)}")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
