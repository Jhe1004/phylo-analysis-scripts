#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/21_annotation/Genoma/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/21_annotation/Genoma/output"
CONDA_ENV_NAME = "gemoma"

TARGET_GENOME_EXTENSION = ".fna"
REFERENCE_GENOME_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/21_annotation/Genoma/input/reference/reference_genome.fasta"
REFERENCE_ANNOTATION_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/21_annotation/Genoma/input/reference/reference_annotation.gff"
THREADS = 70
JAVA_MEMORY = "50G"
GEMOMA_ARGS = [
    "GeMoMaPipeline",
    "GeMoMa.Score=ReAlign",
    "AnnotationFinalizer.r=NO",
    "o=true",
    "pc=true",
]


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
    raise FileNotFoundError(
        f"未找到可执行文件 '{name}'。查找顺序: dependencies/bin -> conda 环境 '{CONDA_ENV_NAME}' -> PATH"
    )


def main() -> None:
    gemoma = find_executable("GeMoMa")
    print(f"使用 GeMoMa: {gemoma}")

    target_genomes = sorted(INPUT_DIR.glob(f"*{TARGET_GENOME_EXTENSION}"))
    if not target_genomes:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 *{TARGET_GENOME_EXTENSION} 目标基因组文件。")

    reference_genome = SCRIPT_DIR / INPUT_DIRECTORY / REFERENCE_GENOME_FILE
    reference_annotation = SCRIPT_DIR / INPUT_DIRECTORY / REFERENCE_ANNOTATION_FILE
    if not reference_genome.exists():
        raise FileNotFoundError(f"未找到参考基因组文件: {reference_genome}")
    if not reference_annotation.exists():
        raise FileNotFoundError(f"未找到参考注释文件: {reference_annotation}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_lines = [
        f"参考基因组: {reference_genome}",
        f"参考注释: {reference_annotation}",
        f"线程数: {THREADS}",
        f"Java 内存: {JAVA_MEMORY}",
        "",
    ]

    for target_genome in target_genomes:
        outdir = OUTPUT_DIR / target_genome.stem
        command = [
            str(gemoma),
            f"-Xmx{JAVA_MEMORY}",
            *GEMOMA_ARGS,
            f"threads={THREADS}",
            f"outdir={outdir}",
            f"t={target_genome}",
            f"a={reference_annotation}",
            f"g={reference_genome}",
        ]
        print("运行命令:", " ".join(map(str, command)))
        subprocess.run(command, check=True)
        summary_lines.append(f"{target_genome.name}\t{outdir.name}")

    (OUTPUT_DIR / "run_summary.tsv").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("GeMoMa 注释完成。")


if __name__ == "__main__":
    main()
