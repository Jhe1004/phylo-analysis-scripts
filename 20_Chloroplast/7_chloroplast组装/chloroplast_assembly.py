#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import os
import shutil
import subprocess
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
HELPER = SCRIPT_DIR / "dependencies" / "scripts" / "run_getorganelle_batch_helper.py"

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

READS_SUBDIRECTORY = "reads"
ASSEMBLY_TYPE = "plastome"
FORWARD_SUFFIX = "_1.clean.fq.gz"
REVERSE_SUFFIX = "_2.clean.fq.gz"
THREADS_PER_JOB = 4
PROCESS_COUNT = 4
# ==============================


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def symlink_or_copy(src: Path, dst: Path) -> None:
    try:
        os.symlink(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def main() -> None:
    input_dir = SCRIPT_DIR / INPUT_DIRECTORY / READS_SUBDIRECTORY
    output_dir = SCRIPT_DIR / OUTPUT_DIRECTORY
    workspace_dir = output_dir / "workspace"
    ensure_directory(output_dir)
    if workspace_dir.exists():
        shutil.rmtree(workspace_dir)
    ensure_directory(workspace_dir)

    for file_path in input_dir.iterdir():
        if file_path.is_file():
            symlink_or_copy(file_path, workspace_dir / file_path.name)

    subprocess.run(
        [
            sys.executable,
            str(HELPER),
            "-t",
            ASSEMBLY_TYPE,
            "-s1",
            FORWARD_SUFFIX,
            "-s2",
            REVERSE_SUFFIX,
            "-T",
            str(THREADS_PER_JOB),
            "-p",
            str(PROCESS_COUNT),
        ],
        cwd=str(workspace_dir),
        check=True,
    )
    print(f"完成。输出目录: {output_dir}")


if __name__ == "__main__":
    main()
