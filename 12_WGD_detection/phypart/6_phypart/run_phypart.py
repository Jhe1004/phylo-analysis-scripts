#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import os
import shutil
import subprocess


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/6_phypart/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/6_phypart/output"

CONDA_ENV_NAME = "trinity_env"
JAVA_EXECUTABLE_NAME = "java"
MAVEN_EXECUTABLE_NAME = "mvn"

SPECIES_TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/6_phypart/input/species_tree.newick"
GENE_TREE_SUBDIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/6_phypart/input/gene_trees"

ANALYSIS_MODE = 2
VERBOSE = True
OUTPUT_PREFIX = "phypart"

AUTO_BUILD_JAR_IF_MISSING = False
PHYPART_JAR = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/6_phypart/dependencies/upstream/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar"
UPSTREAM_SOURCE_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/6_phypart/dependencies/upstream"
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def find_executable(executable_name: str) -> Path | None:
    local_candidate = SCRIPT_DIR / "dependencies" / "bin" / executable_name
    if local_candidate.exists():
        return local_candidate.resolve()
    conda_candidate = Path.home() / "miniconda3" / "envs" / CONDA_ENV_NAME / "bin" / executable_name
    if conda_candidate.exists():
        return conda_candidate.resolve()
    system_candidate = shutil.which(executable_name)
    return Path(system_candidate).resolve() if system_candidate else None


def ensure_phypart_jar() -> Path:
    jar_path = resolve_path(PHYPART_JAR)
    if jar_path.exists():
        return jar_path
    if not AUTO_BUILD_JAR_IF_MISSING:
        raise FileNotFoundError(f"未找到 PhyParts jar: {jar_path}")

    maven_executable = find_executable(MAVEN_EXECUTABLE_NAME)
    if maven_executable is None:
        raise FileNotFoundError(
            f"未找到 {MAVEN_EXECUTABLE_NAME}，无法自动构建 PhyParts jar。请先准备好: {jar_path}"
        )

    upstream_dir = resolve_path(UPSTREAM_SOURCE_DIRECTORY)
    if not (upstream_dir / "pom.xml").exists():
        raise FileNotFoundError(f"未找到 pom.xml，无法自动构建: {upstream_dir / 'pom.xml'}")

    result = subprocess.run(
        [str(maven_executable), "clean", "compile", "assembly:single"],
        cwd=str(upstream_dir),
        capture_output=True,
        text=True,
        check=False,
    )
    build_log = resolve_path(OUTPUT_DIRECTORY) / "phypart_build.log"
    ensure_directory(build_log.parent)
    build_log.write_text((result.stdout or "") + (result.stderr or ""), encoding="utf-8")
    if result.returncode != 0 or not jar_path.exists():
        raise RuntimeError(f"自动构建 PhyParts 失败，请查看日志: {build_log}")
    return jar_path


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    ensure_directory(output_dir)

    java_executable = find_executable(JAVA_EXECUTABLE_NAME)
    if java_executable is None:
        raise FileNotFoundError(f"未找到 {JAVA_EXECUTABLE_NAME}")
    jar_path = ensure_phypart_jar()

    species_tree = input_dir / SPECIES_TREE_FILE
    gene_tree_dir = input_dir / GENE_TREE_SUBDIRECTORY
    if not species_tree.exists():
        raise FileNotFoundError(f"未找到物种树文件: {species_tree}")
    if not gene_tree_dir.exists():
        raise FileNotFoundError(f"未找到基因树目录: {gene_tree_dir}")

    output_prefix_path = output_dir / OUTPUT_PREFIX
    cmd = [
        str(java_executable),
        "-jar",
        str(jar_path),
        "-a",
        str(ANALYSIS_MODE),
        "-d",
        str(gene_tree_dir),
        "-m",
        str(species_tree),
        "-o",
        str(output_prefix_path),
    ]
    if VERBOSE:
        cmd.insert(5, "-v")

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    log_file = output_dir / "phypart.log"
    log_file.write_text((result.stdout or "") + (result.stderr or ""), encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"PhyParts 运行失败，请查看日志: {log_file}")

    summary_file = output_dir / "phypart_summary.txt"
    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write(f"java = {java_executable}\n")
        handle.write(f"jar = {jar_path}\n")
        handle.write(f"species_tree = {species_tree}\n")
        handle.write(f"gene_tree_dir = {gene_tree_dir}\n")
        handle.write(f"analysis_mode = {ANALYSIS_MODE}\n")
        handle.write(f"output_prefix = {output_prefix_path}\n")
        handle.write(f"log_file = {log_file}\n")

    print(f"完成。输出目录: {output_dir}")


if __name__ == "__main__":
    main()
