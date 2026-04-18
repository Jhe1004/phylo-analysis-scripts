#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
HELPER_DIR = SCRIPT_DIR / "dependencies" / "scripts"
sys.path.insert(0, str(HELPER_DIR))

from plot_combined_helper import run_plot_combined
from plot_gmm_helper import run_plot_gmm


# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/absolute_dating/7_visualize_results/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/absolute_dating/7_visualize_results/output"

AGE_SUBDIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/absolute_dating/7_visualize_results/input/ages"
CLADE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/absolute_dating/7_visualize_results/input/clade.md"
TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/absolute_dating/7_visualize_results/input/species_tree.newick"
NAME_MAPPING_FILE = ""

RUN_GMM_PLOTS = True
RUN_COMBINED_PLOTS = True
# ==============================


def main() -> None:
    input_dir = (SCRIPT_DIR / INPUT_DIRECTORY).resolve()
    output_dir = (SCRIPT_DIR / OUTPUT_DIRECTORY).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    ages_dir = input_dir / AGE_SUBDIRECTORY
    clade_file = input_dir / CLADE_FILE
    tree_file = input_dir / TREE_FILE
    mapping_file = (input_dir / NAME_MAPPING_FILE) if NAME_MAPPING_FILE else None

    if RUN_GMM_PLOTS:
        run_plot_gmm(ages_dir=ages_dir, output_dir=output_dir / "gmm_plots")
    if RUN_COMBINED_PLOTS:
        run_plot_combined(
            ages_dir=ages_dir,
            clade_file=clade_file,
            tree_file=tree_file,
            output_dir=output_dir / "combined_output",
            mapping_file=mapping_file,
        )
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
