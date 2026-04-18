#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import shutil
import subprocess
from pathlib import Path


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/9_复杂基因结构建树/1_proteinortho/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/9_复杂基因结构建树/1_proteinortho/output"
CONDA_ENV_NAME = "trinity_env"
INPUT_EXTENSION = ".fasta"
PROJECT_PREFIX = "chloroplast_fragments"


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
    proteinortho = find_executable("proteinortho")
    blastn = find_executable("blastn")
    print(f"使用 proteinortho: {proteinortho}")
    print(f"使用 blastn: {blastn}")

    fasta_files = sorted(INPUT_DIR.glob(f"*{INPUT_EXTENSION}"))
    if not fasta_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 *{INPUT_EXTENSION} 文件。")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    workspace = OUTPUT_DIR / "workspace"
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.mkdir(parents=True, exist_ok=True)

    workspace_files = []
    for fasta_file in fasta_files:
        target = workspace / fasta_file.name
        shutil.copy2(fasta_file, target)
        workspace_files.append(target.name)

    command = [str(proteinortho), "-project", PROJECT_PREFIX, "-p=blastn", *workspace_files]
    subprocess.run(command, cwd=workspace, check=True)

    for item in workspace.iterdir():
        shutil.move(str(item), OUTPUT_DIR / item.name)
    workspace.rmdir()
    print("proteinortho 运行完成。")


if __name__ == "__main__":
    main()
