#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
HELPER_DIR = SCRIPT_DIR / "dependencies" / "scripts"

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/4_高变区/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/4_高变区/output"

ALIGNMENT_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/4_高变区/input/alignment.fasta"
WINDOW_SIZE = 600
STEP_SIZE = 200
OUTPUT_DATA_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/4_高变区/output/pi_window.tsv"
OUTPUT_PLOT_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/4_高变区/output/pi_distribution.png"
PLOT_TITLE = "Sliding Window Analysis of Nucleotide Diversity (π)"
PLOT_COLOR = "#1f77b4"
# ==============================


def main() -> None:
    input_dir = SCRIPT_DIR / INPUT_DIRECTORY
    output_dir = SCRIPT_DIR / OUTPUT_DIRECTORY
    output_dir.mkdir(parents=True, exist_ok=True)
    alignment_file = input_dir / ALIGNMENT_FILE
    data_file = output_dir / OUTPUT_DATA_FILE

    calc = subprocess.run(
        [
            sys.executable,
            str(HELPER_DIR / "calculate_pi_helper.py"),
            str(alignment_file),
            "-w",
            str(WINDOW_SIZE),
            "-s",
            str(STEP_SIZE),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    data_file.write_text(calc.stdout, encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(HELPER_DIR / "plot_pi_helper.py"),
            str(data_file),
            "-o",
            str(output_dir / OUTPUT_PLOT_FILE),
            "--title",
            PLOT_TITLE,
            "--color",
            PLOT_COLOR,
        ],
        check=True,
    )

    print(f"完成。输出目录: {output_dir}")


if __name__ == "__main__":
    main()
