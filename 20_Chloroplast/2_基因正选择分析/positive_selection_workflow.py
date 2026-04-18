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
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/2_基因正选择分析/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/2_基因正选择分析/output"

GENBANK_SUBDIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/2_基因正选择分析/input/gb_files"
TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/2_基因正选择分析/input/species_tree.nwk"
FOREGROUND_SPECIES_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/2_基因正选择分析/input/foreground_species_list.txt"
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

    shutil.copytree(input_dir / GENBANK_SUBDIRECTORY, workspace_dir / "input_gb_files")
    for name, target_name in [(TREE_FILE, "RAxML_bipartitions.result.newick"), (FOREGROUND_SPECIES_FILE, "foreground_species_list.txt")]:
        source = input_dir / name
        if source.exists():
            shutil.copy2(source, workspace_dir / target_name)

    run_helper("extract_genes_helper.py", workspace_dir)
    run_helper("align_and_prep_helper.py", workspace_dir)
    run_helper("prep_paml_files_helper.py", workspace_dir)
    run_helper("run_paml_null_model_helper.py", workspace_dir)
    run_helper("run_and_summarize_paml_helper.py", workspace_dir)
    if RUN_HEATMAP:
        run_helper("plot_heatmap_helper.py", workspace_dir)

    for name in [
        "output_cds_genes",
        "output_pep_genes",
        "aligned_cds",
        "aligned_pep",
        "paml_input",
        "homologs.txt",
        "paml_final_results.csv",
        "kaks_heatmap_final.png",
        "results_summary_english.txt",
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
