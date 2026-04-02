#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量生成 treePL 配置文件。
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import os
import shutil

from ete3 import Tree


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

TREE_SUBDIRECTORY = "trees"
ALIGNMENT_SUBDIRECTORY = "alignments"
CALIBRATED_TREE_FILE = "dna.tree"

TREE_FILE_EXTENSION = ".new"
ALIGNMENT_EXTENSIONS = [".fasta", ".fas", ".fa", ".phy", ".phylip", ".nex", ".nexus"]

TREEPL_THREADS = 12
TREEPL_SMOOTH = 1

PROCESS_COUNT = max(1, min(8, os.cpu_count() or 1))
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def strip_suffix(name: str, suffix: str) -> str:
    return name[: -len(suffix)] if suffix and name.endswith(suffix) else Path(name).stem


def get_alignment_length(alignment_file: Path) -> int:
    if alignment_file.suffix.lower() in {".nex", ".nexus"}:
        with open(alignment_file, "r", encoding="utf-8") as handle:
            for line in handle:
                lower_line = line.strip().lower()
                if "nchar=" in lower_line:
                    chunk = lower_line.split("nchar=", 1)[1]
                    digits = []
                    for char in chunk:
                        if char.isdigit():
                            digits.append(char)
                        elif digits:
                            break
                    if digits:
                        return int("".join(digits))
    sequence_parts: list[str] = []
    with open(alignment_file, "r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith(">"):
                continue
            sequence_parts.append(text)
    return len("".join(sequence_parts[:1])) if sequence_parts else 0


def find_alignment_file(base_name: str, alignment_dir: Path) -> Path | None:
    for extension in ALIGNMENT_EXTENSIONS:
        candidate = alignment_dir / f"{base_name}{extension}"
        if candidate.exists():
            return candidate
    return None


def get_calibrated_age(calibrated_tree: Tree, monophyletic_group: list[str]) -> float:
    if not monophyletic_group:
        return 0.0
    try:
        mrca = calibrated_tree.get_common_ancestor(monophyletic_group)
        any_leaf = mrca.get_leaves()[0]
        return float(mrca.get_distance(any_leaf))
    except Exception:
        return 0.0


def is_monophyletic_in_species_tree(gene_tree: Tree, species_tree: Tree, node) -> bool:
    ingroup = node.get_leaf_names()
    if len(ingroup) <= 1:
        return False
    try:
        species_mrca = species_tree.get_common_ancestor(ingroup)
    except Exception:
        return False
    outgroup = set(gene_tree.get_leaf_names()) - set(ingroup)
    return all(name not in outgroup for name in species_mrca.get_leaf_names())


def generate_calibration_lines(tree_file: Path, calibrated_tree: Tree) -> list[str]:
    gene_tree = Tree(str(tree_file), format=1)
    calibrations: list[str] = []

    all_leaves = gene_tree.get_leaf_names()
    root_age = get_calibrated_age(calibrated_tree, all_leaves)
    if root_age > 0:
        calibrations.append(f"mrca = root_node {' '.join(all_leaves)}\n")
        calibrations.append(f"fixage = root_node {root_age}\n")

    node_counter = 1
    for node in gene_tree.traverse():
        if node.is_leaf():
            continue
        if not is_monophyletic_in_species_tree(gene_tree, calibrated_tree, node):
            continue
        node_leaves = node.get_leaf_names()
        node_age = get_calibrated_age(calibrated_tree, node_leaves)
        if node_age <= 0:
            continue
        node_name = f"node_{node_counter}"
        calibrations.append(f"mrca = {node_name} {' '.join(node_leaves)}\n")
        calibrations.append(f"min = {node_name} {node_age * 0.95}\n")
        calibrations.append(f"max = {node_name} {node_age * 1.05}\n")
        node_counter += 1
    return calibrations


def build_conf_text(
    tree_file: Path,
    alignment_file: Path | None,
    calibrated_tree: Tree,
    output_base_dir: Path,
) -> str:
    config_dir = output_base_dir / "configs"
    work_tree_dir = output_base_dir / "trees"
    dated_tree_dir = output_base_dir / "dated_trees"
    log_dir = output_base_dir / "logs"

    relative_tree_path = os.path.relpath(work_tree_dir / tree_file.name, config_dir)
    relative_outfile = os.path.relpath(dated_tree_dir / f"{tree_file.name}.newtree", config_dir)
    relative_cvoutfile = os.path.relpath(log_dir / f"{tree_file.name}.cvout", config_dir)

    lines = [
        f"// Configuration file for {tree_file.name}\n\n",
        "[Input files containing the ML trees]\n",
        f"treefile = {relative_tree_path}\n\n",
        "[General commands]\n",
        f"nthreads = {TREEPL_THREADS}\n",
        f"smooth = {TREEPL_SMOOTH}\n",
        "thorough\n",
        "log_pen\n",
    ]

    if alignment_file is not None:
        lines.append(f"numsites = {get_alignment_length(alignment_file)}\n\n")
    else:
        lines.append("// numsites = 0  请手动补充序列长度\n\n")

    lines.append("[Calibrations]\n")
    calibration_lines = generate_calibration_lines(tree_file, calibrated_tree)
    if calibration_lines:
        lines.extend(calibration_lines)
    else:
        lines.append("// No valid calibration points found.\n")
    lines.append("\n")

    lines.extend(
        [
            "[Optimisation parameters]\n",
            "opt = 3\n",
            "moredetail\n",
            "optad = 3\n",
            "moredetailad\n",
            "optcvad = 5\n",
            "randomcv\n",
            "cviter = 5\n",
            "cvsimaniter = 100000\n",
            "cvstart = 10000\n",
            "cvstop = 0.00001\n",
            "cvmultstep = 0.1\n",
            f"cvoutfile = {relative_cvoutfile}\n",
            f"outfile = {relative_outfile}\n",
        ]
    )
    return "".join(lines)


def generate_single_config(task: tuple[Path, Path, Path, str]) -> tuple[str, str]:
    tree_file, alignment_dir, output_dir, calibrated_tree_file = task
    calibrated_tree = Tree(str(calibrated_tree_file), format=1)
    base_name = strip_suffix(tree_file.name, TREE_FILE_EXTENSION)
    alignment_file = find_alignment_file(base_name, alignment_dir)
    config_text = build_conf_text(tree_file, alignment_file, calibrated_tree, output_dir)
    config_path = output_dir / "configs" / f"{base_name}.conf"
    with open(config_path, "w", encoding="utf-8") as handle:
        handle.write(config_text)
    return tree_file.name, config_path.name


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    tree_input_dir = input_dir / TREE_SUBDIRECTORY
    alignment_input_dir = input_dir / ALIGNMENT_SUBDIRECTORY
    calibrated_tree_file = input_dir / CALIBRATED_TREE_FILE

    config_output_dir = output_dir / "configs"
    copied_tree_dir = output_dir / "trees"
    dated_tree_dir = output_dir / "dated_trees"
    log_dir = output_dir / "logs"

    for path in [config_output_dir, copied_tree_dir, dated_tree_dir, log_dir]:
        ensure_directory(path)

    tree_files = sorted(tree_input_dir.glob(f"*{TREE_FILE_EXTENSION}"))
    if not tree_files:
        raise FileNotFoundError(f"未找到树文件: {tree_input_dir}")
    if not calibrated_tree_file.exists():
        raise FileNotFoundError(f"未找到校准物种树: {calibrated_tree_file}")

    for tree_file in tree_files:
        shutil.copy2(tree_file, copied_tree_dir / tree_file.name)

    tasks = [
        (copied_tree_dir / tree_file.name, alignment_input_dir, output_dir, calibrated_tree_file)
        for tree_file in tree_files
    ]

    if PROCESS_COUNT > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=PROCESS_COUNT) as executor:
            results = list(executor.map(generate_single_config, tasks))
    else:
        results = [generate_single_config(task) for task in tasks]

    print(f"完成。共生成 {len(results)} 个 treePL 配置文件。")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
