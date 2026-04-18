#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
按顺序完成 DensiTree 绘图、节点支持统计和支持率饼图绘制。
"""

from __future__ import annotations

from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
HELPER_DIR = SCRIPT_DIR / "dependencies" / "scripts"
sys.path.insert(0, str(HELPER_DIR))

from count_species_tree_node_support_helper import run_count_species_tree_node_support
from densitree_plot_helper import run_densitree_plot
from draw_species_tree_support_pies_helper import run_draw_species_tree_support_pies


# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/3_draw_densitree/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/3_draw_densitree/output"

GENE_TREES_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/3_draw_densitree/input/result.trees"
SPECIES_TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/3_draw_densitree/input/species_tree.newick"
TIP_ORDER_FILE = ""
NAME_MAPPING_FILE = ""

RUN_DENSITREE_PLOT = True
RUN_SUPPORT_COUNT = True
RUN_SUPPORT_PIE_PLOT = True

DENSITREE_OUTPUT_FILE = "densi_tree_plot.svg"
SUPPORT_CSV_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/3_draw_densitree/output/species_tree_node_support.csv"
SUPPORT_PIE_OUTPUT_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/3_draw_densitree/output/species_tree_support_pies.svg"

PROCESS_COUNT = 4
IMAGE_WIDTH = 1000
IMAGE_HEIGHT = 1200
PIE_MAX_RADIUS = 11.0
PIE_SCALE = 0.30
# ==============================


def build_path(relative_directory: str, file_name: str) -> Path:
    return (SCRIPT_DIR / relative_directory / file_name).resolve()


def main() -> None:
    input_dir = (SCRIPT_DIR / INPUT_DIRECTORY).resolve()
    output_dir = (SCRIPT_DIR / OUTPUT_DIRECTORY).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    species_tree_path = build_path(INPUT_DIRECTORY, SPECIES_TREE_FILE)
    gene_trees_path = build_path(INPUT_DIRECTORY, GENE_TREES_FILE)
    tip_order_path = build_path(INPUT_DIRECTORY, TIP_ORDER_FILE) if TIP_ORDER_FILE else None
    name_mapping_path = build_path(INPUT_DIRECTORY, NAME_MAPPING_FILE) if NAME_MAPPING_FILE else None
    support_csv_path = build_path(OUTPUT_DIRECTORY, SUPPORT_CSV_FILE)

    if RUN_DENSITREE_PLOT:
        run_densitree_plot(
            species_tree_path=species_tree_path,
            gene_trees_path=gene_trees_path,
            output_path=build_path(OUTPUT_DIRECTORY, DENSITREE_OUTPUT_FILE),
            tip_order_path=tip_order_path,
            name_mapping_path=name_mapping_path,
            width=IMAGE_WIDTH,
            height=IMAGE_HEIGHT,
        )

    if RUN_SUPPORT_COUNT:
        run_count_species_tree_node_support(
            species_tree_path=species_tree_path,
            gene_trees_path=gene_trees_path,
            output_path=support_csv_path,
            name_mapping_path=name_mapping_path,
            process_count=PROCESS_COUNT,
        )

    if RUN_SUPPORT_PIE_PLOT:
        run_draw_species_tree_support_pies(
            species_tree_path=species_tree_path,
            gene_trees_path=gene_trees_path,
            support_csv_path=support_csv_path,
            output_path=build_path(OUTPUT_DIRECTORY, SUPPORT_PIE_OUTPUT_FILE),
            tip_order_path=tip_order_path,
            name_mapping_path=name_mapping_path,
            width=IMAGE_WIDTH,
            height=IMAGE_HEIGHT,
            pie_max_radius=PIE_MAX_RADIUS,
            pie_scale=PIE_SCALE,
        )

    print(f"完成。输出目录: {output_dir}")


if __name__ == "__main__":
    main()
