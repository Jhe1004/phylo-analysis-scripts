#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
from pathlib import Path

import toyplot
import toyplot.svg
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


def load_tip_order(tip_order_file: Path | None, name_mapping: dict[str, str]) -> list[str] | None:
    if tip_order_file is None or not tip_order_file.exists():
        return None
    with open(tip_order_file, "r", encoding="utf-8") as handle:
        raw_order = [line.strip() for line in handle if line.strip()]
    return [name_mapping.get(name, name) for name in raw_order[::-1]]


def validate_tip_order(tree, tip_order: list[str] | None) -> None:
    if tip_order is None:
        return
    tree_tips = set(tree.get_tip_labels())
    missing = [name for name in tip_order if name not in tree_tips]
    if missing:
        raise ValueError(f"顺序文件中的名称与树不匹配，前 10 个缺失名称: {missing[:10]}")


def scale_tree_to_height(tree, target_height: float = 1.0):
    current_height = tree.treenode.height
    if current_height == 0:
        return tree
    scaling_factor = target_height / current_height
    for node in tree.treenode.traverse():
        node._dist *= scaling_factor
    return tree


def run_densitree_plot(
    species_tree_path: Path,
    gene_trees_path: Path,
    output_path: Path,
    tip_order_path: Path | None = None,
    name_mapping_path: Path | None = None,
    width: int = 1000,
    height: int = 1200,
) -> None:
    print(f"正在加载物种树: {species_tree_path}")
    species_tree = toytree.tree(str(species_tree_path))
    print(f"正在加载基因树集: {gene_trees_path}")
    gene_trees = toytree.mtree(str(gene_trees_path))

    name_mapping = load_name_mapping(name_mapping_path)
    species_tree = rename_tree_tips(species_tree, name_mapping)
    gene_trees.treelist = [rename_tree_tips(tree, name_mapping) for tree in gene_trees]

    tip_order = load_tip_order(tip_order_path, name_mapping)
    validate_tip_order(species_tree, tip_order)

    species_tree = scale_tree_to_height(species_tree, 1.0)
    gene_trees.treelist = [scale_tree_to_height(tree, 1.0) for tree in gene_trees]

    canvas = toyplot.Canvas(width=width, height=height)
    axes = canvas.cartesian(xlabel="Relative time (root to tip)")

    gene_trees.draw_cloud_tree(
        axes=axes,
        fixed_order=tip_order,
        edge_style={
            "stroke": "#4C72B0",
            "stroke-opacity": 0.10,
            "stroke-width": 1.0,
        },
    )
    species_tree.draw(
        axes=axes,
        fixed_order=tip_order,
        edge_type="c",
        edge_style={
            "stroke": "black",
            "stroke-width": 1.0,
        },
        tip_labels_align=True,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    toyplot.svg.render(canvas, str(output_path))
    print(f"已写出 DensiTree 图: {output_path}")
