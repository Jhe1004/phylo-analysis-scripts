#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import os
import re


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/5_rename/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/5_rename/output"

INPUT_EXTENSION = ".tree"
OUTPUT_SUFFIX = "_re.tree"
PROCESS_COUNT = max(1, min(8, os.cpu_count() or 1))
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def clean_tree_labels(text: str) -> str:
    return re.sub(r"\+\+[^,:;()]+", "", text)


def process_tree(task: tuple[Path, Path]) -> str:
    input_file, output_file = task
    output_file.write_text(clean_tree_labels(input_file.read_text(encoding="utf-8")), encoding="utf-8")
    return input_file.name


def process_tree_task(task: tuple[Path, Path]) -> str:
    return process_tree(task)


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    ensure_directory(output_dir)

    tree_files = sorted(input_dir.glob(f"*{INPUT_EXTENSION}"))
    if not tree_files:
        raise FileNotFoundError(f"未在 {input_dir} 中找到 {INPUT_EXTENSION} 文件")

    tasks = [(tree_file, output_dir / f"{tree_file.stem}{OUTPUT_SUFFIX}") for tree_file in tree_files]
    if PROCESS_COUNT > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=PROCESS_COUNT) as executor:
            results = list(executor.map(process_tree_task, tasks))
    else:
        results = [process_tree_task(task) for task in tasks]

    summary_file = output_dir / "rename_summary.tsv"
    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write("tree\tstatus\n")
        for name in results:
            handle.write(f"{name}\tsuccess\n")

    print(f"完成树数: {len(results)}")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
