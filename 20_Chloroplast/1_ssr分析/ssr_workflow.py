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
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/1_ssr分析/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/1_ssr分析/output"

GENBANK_SUBDIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/1_ssr分析/input/gb_files"
SPECIES_TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/1_ssr分析/input/species_tree.nwk"
MISA_SCRIPT_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/1_ssr分析/dependencies/tools/misa.py"
MISA_CONFIG_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/1_ssr分析/dependencies/tools/misa.ini"
RUN_SUMMARY_PLOT = True
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

    shutil.copytree(input_dir / GENBANK_SUBDIRECTORY, workspace_dir / "input_gb_files")
    species_tree_source = input_dir / SPECIES_TREE_FILE
    if species_tree_source.exists():
        shutil.copy2(species_tree_source, workspace_dir / species_tree_source.name)
    for tool_name in [MISA_SCRIPT_FILE, MISA_CONFIG_FILE]:
        tool_source = SCRIPT_DIR / "dependencies" / "tools" / tool_name
        if tool_source.exists():
            shutil.copy2(tool_source, workspace_dir / tool_source.name)

    run_helper("gb_to_fasta_helper.py", workspace_dir)
    run_helper("find_ssrs_helper.py", workspace_dir)
    if RUN_SUMMARY_PLOT:
        run_helper("summarize_ssrs_helper.py", workspace_dir)

    for name in [
        "output_genome_fasta",
        "output_ssr_results_misa",
        "ssr_summary_data.csv",
        "ssr_results_report.txt",
        "ssr_phylo_stacked_barplot_corrected.png",
    ]:
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
