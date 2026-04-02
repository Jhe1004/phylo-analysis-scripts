#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
CONDA_ENV_NAME = "trinity_env"
PROCESS_COUNT = 8
THREADS_PER_PROCESS = 1
GROUP_DIRECTORY_GLOB = "*_groups"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY
HELPER_DIR = SCRIPT_DIR / "dependencies" / "scripts"
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
    mafft = find_executable("mafft")
    print(f"使用 mafft: {mafft}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    workspace = OUTPUT_DIR / "workspace"
    if workspace.exists():
        shutil.rmtree(workspace)
    shutil.copytree(INPUT_DIR, workspace)

    env = dict(**os.environ)
    env["PATH"] = f"{mafft.parent}:{env.get('PATH', '')}"
    env["NUM_CPU"] = str(PROCESS_COUNT)
    env["THREAD_PER_JOB"] = str(THREADS_PER_PROCESS)
    env["SUBDIR_PATTERN"] = GROUP_DIRECTORY_GLOB

    subprocess.run([sys.executable, str(HELPER_DIR / "mafft_batch_helper.py")], cwd=workspace, env=env, check=True)
    subprocess.run([sys.executable, str(HELPER_DIR / "concat_batch_helper.py")], cwd=workspace, env=env, check=True)
    print(f"MAFFT 和拼接完成，结果在 {workspace}")


if __name__ == "__main__":
    import os

    main()
