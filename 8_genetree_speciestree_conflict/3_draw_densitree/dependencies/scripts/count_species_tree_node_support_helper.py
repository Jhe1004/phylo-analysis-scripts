#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import csv
from pathlib import Path

import toytree


def load_name_mapping(mapping_file: Path | None) -> dict[str, str]:
    if mapping_file is None or not mapping_file.exists():
        return {}
    mapping = {}
    with open(mapping_file, "r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            current_name = row["current_sample_name"].strip()
            standard_name = row["standard_latin_name"].strip()
            if current_name and standard_name:
                mapping[current_name] = standard_name
    return mapping


def rename_tree_tips(tree, name_mapping: dict[str, str]):
    if not name_mapping:
        return tree
    for node in tree.treenode.traverse():
        if node.is_leaf() and node.name in name_mapping:
            node.name = name_mapping[node.name]
    return tree


def node_clade(node) -> frozenset[str]:
    return frozenset(leaf.name for leaf in node.iter_leaves())


def collect_internal_clades(tree) -> set[frozenset[str]]:
    return {
        node_clade(node)
        for node in tree.treenode.traverse()
        if not node.is_leaf()
    }


def compute_gene_clades(task: tuple[str, dict[str, str]]) -> tuple[set[frozenset[str]], set[str]]:
    newick, name_mapping = task
    tree = toytree.tree(newick)
    tree = rename_tree_tips(tree, name_mapping)
    return collect_internal_clades(tree), set(tree.get_tip_labels())


def run_count_species_tree_node_support(
    species_tree_path: Path,
    gene_trees_path: Path,
    output_path: Path,
    name_mapping_path: Path | None = None,
    process_count: int = 1,
) -> None:
    print(f"正在加载物种树: {species_tree_path}")
    species_tree = toytree.tree(str(species_tree_path))
    name_mapping = load_name_mapping(name_mapping_path)
    species_tree = rename_tree_tips(species_tree, name_mapping)

    with open(gene_trees_path, "r", encoding="utf-8") as handle:
        gene_tree_strings = [line.strip() for line in handle if line.strip()]
    if not gene_tree_strings:
        raise ValueError(f"基因树文件为空: {gene_trees_path}")

    tasks = [(newick, name_mapping) for newick in gene_tree_strings]
    if process_count > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=process_count) as executor:
            gene_results = list(executor.map(compute_gene_clades, tasks))
    else:
        gene_results = [compute_gene_clades(task) for task in tasks]

    species_tips = set(species_tree.get_tip_labels())
    for index, (_, gene_tips) in enumerate(gene_results, start=1):
        if gene_tips != species_tips:
            missing = sorted(species_tips - gene_tips)
            extra = sorted(gene_tips - species_tips)
            raise ValueError(
                f"第 {index} 棵基因树的 tip 集与物种树不一致。"
                f" missing={missing[:10]} extra={extra[:10]}"
            )

    gene_clade_sets = [item[0] for item in gene_results]
    internal_nodes = [node for node in species_tree.treenode.traverse() if not node.is_leaf()]
    total_gene_trees = len(gene_clade_sets)

    rows = []
    for order, node in enumerate(internal_nodes, start=1):
        clade = node_clade(node)
        support_count = sum(1 for clade_set in gene_clade_sets if clade in clade_set)
        support_percent = 100.0 * support_count / total_gene_trees if total_gene_trees else 0.0
        rows.append(
            {
                "node_order": order,
                "node_idx": getattr(node, "idx", getattr(node, "_idx", order)),
                "clade_size": len(clade),
                "supporting_gene_trees": support_count,
                "support_percent": f"{support_percent:.2f}",
                "clade_members": ";".join(sorted(clade)),
                "is_root": "yes" if node.is_root() else "no",
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "node_order",
                "node_idx",
                "clade_size",
                "supporting_gene_trees",
                "support_percent",
                "clade_members",
                "is_root",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"已写出节点支持统计: {output_path}")
