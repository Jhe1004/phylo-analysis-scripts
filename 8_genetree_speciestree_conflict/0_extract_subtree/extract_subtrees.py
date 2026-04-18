#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量提取基因树中的目标子树，并统一重新置根。
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import os

from ete3 import Tree


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/0_extract_subtree/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/0_extract_subtree/output"

TREE_SUBDIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/0_extract_subtree/input/trees"
SUBTREE_SPECIES_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_genetree_speciestree_conflict/0_extract_subtree/input/subtree_species.txt"

TREE_FILE_PATTERNS = ["RAxML*"]
OUTPUT_SUFFIX = ".new"

OUTGROUP_SPECIES = ""
INCLUDE_OUTGROUP_IN_SUBTREE = True
REQUIRE_ALL_TARGET_SPECIES = True

PROCESS_COUNT = max(1, min(8, os.cpu_count() or 1))
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_species_list(species_file: Path) -> list[str]:
    species = []
    with open(species_file, "r", encoding="utf-8") as handle:
        for line in handle:
            name = line.strip()
            if name:
                species.append(name)
    if not species:
        raise ValueError(f"目标物种列表为空: {species_file}")
    return species


def discover_tree_files(tree_dir: Path) -> list[Path]:
    files: list[Path] = []
    seen: set[Path] = set()
    for pattern in TREE_FILE_PATTERNS:
        for path in sorted(tree_dir.glob(pattern)):
            if path.is_file() and path not in seen:
                files.append(path)
                seen.add(path)
    return files


def process_tree_file(
    tree_path: Path,
    output_path: Path,
    subtree_species_list: list[str],
    outgroup_species: str,
    require_all_species: bool,
) -> tuple[str, int, int]:
    processed_count = 0
    kept_count = 0
    output_lines: list[str] = []

    with open(tree_path, "r", encoding="utf-8") as handle:
        for line in handle:
            newick = line.strip()
            if not newick:
                continue
            processed_count += 1
            try:
                tree = Tree(newick)
                if require_all_species:
                    missing = [name for name in subtree_species_list if name not in tree.get_leaf_names()]
                    if missing:
                        continue
                tree.prune(subtree_species_list)
                if outgroup_species:
                    if outgroup_species not in tree.get_leaf_names():
                        continue
                    tree.set_outgroup(outgroup_species)
                output_lines.append(tree.write() + "\n")
                kept_count += 1
            except Exception:
                continue

    with open(output_path, "w", encoding="utf-8") as handle:
        handle.writelines(output_lines)

    return tree_path.name, processed_count, kept_count


def process_tree_file_task(task: tuple[Path, Path, list[str], str, bool]) -> tuple[str, int, int]:
    return process_tree_file(*task)


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    tree_dir = input_dir / TREE_SUBDIRECTORY
    species_file = input_dir / SUBTREE_SPECIES_FILE
    output_tree_dir = output_dir / TREE_SUBDIRECTORY
    summary_file = output_dir / "subtree_extraction_summary.tsv"

    ensure_directory(output_tree_dir)

    subtree_species = load_species_list(species_file)
    if OUTGROUP_SPECIES and INCLUDE_OUTGROUP_IN_SUBTREE and OUTGROUP_SPECIES not in subtree_species:
        subtree_species.append(OUTGROUP_SPECIES)

    tree_files = discover_tree_files(tree_dir)
    if not tree_files:
        raise FileNotFoundError(f"未找到待处理树文件: {tree_dir}")

    print(f"找到 {len(tree_files)} 个树文件，开始提取子树。")

    tasks = []
    for tree_path in tree_files:
        output_path = output_tree_dir / f"{tree_path.name}{OUTPUT_SUFFIX}"
        tasks.append((tree_path, output_path, subtree_species, OUTGROUP_SPECIES, REQUIRE_ALL_TARGET_SPECIES))

    if PROCESS_COUNT > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=PROCESS_COUNT) as executor:
            results = list(executor.map(process_tree_file_task, tasks))
    else:
        results = [process_tree_file(*task) for task in tasks]

    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write("tree_file\tprocessed_trees\tkept_trees\n")
        for tree_name, processed_count, kept_count in results:
            handle.write(f"{tree_name}\t{processed_count}\t{kept_count}\n")

    total_processed = sum(item[1] for item in results)
    total_kept = sum(item[2] for item in results)
    print(f"完成。总树数: {total_processed}，成功保留: {total_kept}")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
