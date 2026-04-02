#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
CONDA_ENV_NAME = "trinity_env"
ALIGNMENT_FILE = "my_alignment.fasta"
PARTITION_FILE = "partitions.txt"
IQTREE_EXECUTABLE = "iqtree"
IQTREE_ARGS = ["-m", "MFP", "-bb", "1000", "-nt", "AUTO"]


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY
DEPENDENCIES_DIR = SCRIPT_DIR / "dependencies"


def resolve_conda_executable(env_name: str, executable_name: str) -> Path | None:
    for prefix in (
        Path.home() / "miniconda3" / "envs" / env_name / "bin",
        Path.home() / "mambaforge" / "envs" / env_name / "bin",
        Path.home() / "anaconda3" / "envs" / env_name / "bin",
    ):
        candidate = prefix / executable_name
        if candidate.exists():
            return candidate
    return None


def find_executable(name: str) -> Path:
    local_candidate = DEPENDENCIES_DIR / "bin" / name
    if local_candidate.exists():
        return local_candidate
    conda_candidate = resolve_conda_executable(CONDA_ENV_NAME, name)
    if conda_candidate is not None:
        return conda_candidate
    system_candidate = shutil.which(name)
    if system_candidate:
        return Path(system_candidate)
    raise FileNotFoundError(f"未找到 {name}")


def main() -> None:
    iqtree = find_executable(IQTREE_EXECUTABLE)
    alignment_path = INPUT_DIR / ALIGNMENT_FILE
    partition_path = INPUT_DIR / PARTITION_FILE
    if not alignment_path.exists():
        raise FileNotFoundError(f"未找到比对文件: {alignment_path}")
    if not partition_path.exists():
        raise FileNotFoundError(f"未找到分区文件: {partition_path}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    command = [
        str(iqtree),
        "-s",
        str(alignment_path),
        "-spp",
        str(partition_path),
        *IQTREE_ARGS,
        "-pre",
        str(OUTPUT_DIR / alignment_path.stem),
    ]
    subprocess.run(command, check=True)
    print("IQ-TREE 运行完成。")


if __name__ == "__main__":
    main()
