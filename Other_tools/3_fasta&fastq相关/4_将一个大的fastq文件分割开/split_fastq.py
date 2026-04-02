#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import shutil
import subprocess
from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
CONDA_ENV_NAME = "trinity_env"
SUFFIX1 = "_1.clean.fq.gz"
SUFFIX2 = "_2.clean.fq.gz"
SEQKIT_THREADS = 4
PARTS = 2


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
    raise FileNotFoundError(f"未找到可执行文件: {name}")


def main() -> None:
    seqkit = find_executable("seqkit")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    file_list = sorted(INPUT_DIR.glob(f"*{SUFFIX1}"))
    if not file_list:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 *{SUFFIX1} 文件。")
    for r1_file in file_list:
        prefix = r1_file.name[: -len(SUFFIX1)]
        r2_file = INPUT_DIR / f"{prefix}{SUFFIX2}"
        if not r2_file.exists():
            continue
        command = [
            str(seqkit),
            "split2",
            "-1",
            str(r1_file),
            "-2",
            str(r2_file),
            "-p",
            str(PARTS),
            "-j",
            str(SEQKIT_THREADS),
            "-O",
            str(OUTPUT_DIR),
        ]
        subprocess.run(command, check=True)
    print("FASTQ 分割完成。")


if __name__ == "__main__":
    main()
