#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import random
import shutil
import subprocess
from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
CONDA_ENV_NAME = "trinity_env"
FASTA_ALIGNMENT_FILE = ""
PHYLIP_ALIGNMENT_FILE = "cds_groups.phy"
CONFIG_FILE = "config.nex"
PARTITION_FILE = "partitions.txt"
USE_PARTITION = False
MPI_PROCESS_COUNT = 40
RUN_NAME = "cds"
EXABAYES_EXECUTABLE = "exabayes"
MPIRUN_EXECUTABLE = "mpirun"


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


def ensure_phylip_alignment(workspace: Path) -> Path:
    phylip_path = INPUT_DIR / PHYLIP_ALIGNMENT_FILE
    if phylip_path.exists():
        shutil.copy2(phylip_path, workspace / phylip_path.name)
        return workspace / phylip_path.name
    if not FASTA_ALIGNMENT_FILE:
        raise FileNotFoundError("未找到 PHYLIP 比对文件，且 FASTA_ALIGNMENT_FILE 为空。")
    fasta_path = INPUT_DIR / FASTA_ALIGNMENT_FILE
    if not fasta_path.exists():
        raise FileNotFoundError(f"未找到 FASTA 比对文件: {fasta_path}")
    helper = DEPENDENCIES_DIR / "scripts" / "fasta_to_phylip.py"
    output_path = workspace / Path(PHYLIP_ALIGNMENT_FILE).name
    subprocess.run(["python", str(helper), str(fasta_path), str(output_path)], check=True)
    return output_path


def main() -> None:
    mpirun = find_executable(MPIRUN_EXECUTABLE)
    exabayes = find_executable(EXABAYES_EXECUTABLE)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    workspace = OUTPUT_DIR / "workspace"
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.mkdir(parents=True, exist_ok=True)

    phylip_alignment = ensure_phylip_alignment(workspace)
    config_path = INPUT_DIR / CONFIG_FILE
    if not config_path.exists():
        raise FileNotFoundError(f"未找到配置文件: {config_path}")
    shutil.copy2(config_path, workspace / config_path.name)

    command = [
        str(mpirun),
        "-np",
        str(MPI_PROCESS_COUNT),
        str(exabayes),
        "-f",
        phylip_alignment.name,
        "-m",
        "DNA",
        "-n",
        RUN_NAME,
        "-s",
        str(random.randint(0, 32767)),
        "-c",
        config_path.name,
    ]

    partition_path = INPUT_DIR / PARTITION_FILE
    if USE_PARTITION:
        if not partition_path.exists():
            raise FileNotFoundError(f"未找到分区文件: {partition_path}")
        shutil.copy2(partition_path, workspace / partition_path.name)
        command.extend(["-q", partition_path.name])

    subprocess.run(command, cwd=workspace, check=True)
    print("ExaBayes 运行完成。")


if __name__ == "__main__":
    main()
