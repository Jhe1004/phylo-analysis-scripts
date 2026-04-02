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
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

GENBANK_SUBDIRECTORY = "gb_files"
SPECIES_TREE_FILE = "species_tree.nwk"
RUN_HEATMAP = True
# ==============================


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def run_helper(helper_name: str, cwd: Path) -> None:
    subprocess.run([sys.executable, str(HELPER_DIR / helper_name)], cwd=str(cwd), check=True)


def main() -> None:
    input_dir = SCRIPT_DIR / INPUT_DIRECTORY
    output_dir = SCRIPT_DIR / OUTPUT_DIRECTORY
    workspace_dir = output_dir / "workspace"
    ensure_directory(output_dir)
    if workspace_dir.exists():
        shutil.rmtree(workspace_dir)
    ensure_directory(workspace_dir)

    source_gb_dir = input_dir / GENBANK_SUBDIRECTORY
    target_gb_dir = workspace_dir / "input_gb_files"
    shutil.copytree(source_gb_dir, target_gb_dir)

    tree_file = input_dir / SPECIES_TREE_FILE
    if RUN_HEATMAP and tree_file.exists():
        shutil.copy2(tree_file, workspace_dir / "species_tree.nwk")

    run_helper("extract_cds_helper.py", workspace_dir)
    run_helper("calculate_rscu_helper.py", workspace_dir)
    if RUN_HEATMAP:
        run_helper("plot_heatmap_helper.py", workspace_dir)

    for name in ["output_cds_fasta", "rscu_matrix.csv", "codon_usage_report.txt", "final_heatmap_customized.png"]:
        path = workspace_dir / name
        if path.exists():
            target = output_dir / name
            if path.is_dir():
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(path, target)
            else:
                shutil.copy2(path, target)

    print(f"完成。输出目录: {output_dir}")


if __name__ == "__main__":
    main()
