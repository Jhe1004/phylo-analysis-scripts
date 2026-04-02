#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


# =============================
# 用户配置区
# =============================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
CONDA_ENV_NAME = "reseq"

REFERENCE_GENOME = "reference/ref.fasta"
FASTQ_SUBDIRECTORY = "fastq"

STEPS = {
    "step1_index": False,
    "step2_mapping": True,
    "step3_calling": True,
    "step4_consensus": True,
    "step5_qc": True,
    "step6_combine": True,
}

STEP2_THREADS_PER_SAMPLE = 20
STEP2_PARALLEL_SAMPLES = 10
STEP2_PLUS_TAG = "_1.clean.fq.gz"
STEP2_MINUS_TAG = "_2.clean.fq.gz"

STEP3_PARALLEL_SAMPLES = 20
STEP3_GATK_THREADS = 6
STEP3_JOINT_CALLING = False
STEP3_SINGLE_VCF = False

STEP4_PARALLEL_SAMPLES = 10
STEP4_MIN_GQ = 20
STEP4_MIN_DP = 3
STEP4_USE_IUPAC = True


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY
DEPENDENCIES_DIR = SCRIPT_DIR / "dependencies"
STEP_SCRIPT_DIR = DEPENDENCIES_DIR / "scripts"


def resolve_conda_bin(env_name: str) -> Path | None:
    for prefix in (
        Path.home() / "miniconda3" / "envs" / env_name / "bin",
        Path.home() / "mambaforge" / "envs" / env_name / "bin",
        Path.home() / "anaconda3" / "envs" / env_name / "bin",
        Path("/home/hejian/anaconda3/envs") / env_name / "bin",
    ):
        if prefix.exists():
            return prefix
    return None


def build_runtime_env() -> dict[str, str]:
    env = os.environ.copy()
    path_parts = []
    local_bin = DEPENDENCIES_DIR / "bin"
    if local_bin.exists():
        path_parts.append(str(local_bin))
    conda_bin = resolve_conda_bin(CONDA_ENV_NAME)
    if conda_bin is not None:
        path_parts.append(str(conda_bin))
    path_parts.append(env.get("PATH", ""))
    env["PATH"] = ":".join(part for part in path_parts if part)
    return env


def run_python_script(script_name: str, args: list[str], cwd: Path, env: dict[str, str]) -> None:
    script_path = STEP_SCRIPT_DIR / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"未找到依赖脚本: {script_path}")
    command = [sys.executable, str(script_path), *args]
    print(f"\n[执行] {' '.join(command)}")
    subprocess.run(command, cwd=str(cwd), env=env, check=True)


def run_shell_command(command: str, cwd: Path, env: dict[str, str]) -> None:
    print(f"\n[执行] {command}")
    subprocess.run(command, cwd=str(cwd), env=env, shell=True, check=True)


def main() -> None:
    reference_genome = INPUT_DIR / REFERENCE_GENOME
    fastq_dir = INPUT_DIR / FASTQ_SUBDIRECTORY
    if not reference_genome.exists():
        raise FileNotFoundError(f"未找到参考基因组文件: {reference_genome}")
    if not fastq_dir.exists():
        raise FileNotFoundError(f"未找到 FASTQ 输入目录: {fastq_dir}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    bam_dir = OUTPUT_DIR / "bam_output"
    vcf_dir = OUTPUT_DIR / "vcf_output"
    consensus_dir = OUTPUT_DIR / "consensus_fasta_output"
    for directory in [bam_dir, vcf_dir, consensus_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    env = build_runtime_env()
    start_time = time.time()
    print("=" * 60)
    print(f"重测序流程启动时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    if STEPS["step1_index"]:
        print("\n>>> [Step 1/6] Indexing Reference...")
        run_python_script(
            "1_index_reference.py",
            ["-r", str(reference_genome)],
            cwd=OUTPUT_DIR,
            env=env,
        )

    if STEPS["step2_mapping"]:
        print("\n>>> [Step 2/6] BWA Mapping...")
        run_python_script(
            "2_bwa_map.py",
            [
                "-r",
                str(reference_genome),
                "-f",
                str(fastq_dir),
                "-o",
                str(bam_dir),
                "-p",
                STEP2_PLUS_TAG,
                "-m",
                STEP2_MINUS_TAG,
                "-t",
                str(STEP2_THREADS_PER_SAMPLE),
                "-j",
                str(STEP2_PARALLEL_SAMPLES),
            ],
            cwd=OUTPUT_DIR,
            env=env,
        )

    if STEPS["step3_calling"]:
        print("\n>>> [Step 3/6] GATK Variant Calling...")
        args = [
            "-r",
            str(reference_genome),
            "-b",
            str(bam_dir),
            "-o",
            str(vcf_dir),
            "-t",
            str(STEP3_PARALLEL_SAMPLES),
            "--gatk-threads",
            str(STEP3_GATK_THREADS),
        ]
        if STEP3_JOINT_CALLING:
            args.append("--joint-calling")
        if STEP3_SINGLE_VCF:
            args.append("--single-vcf")
        run_python_script("3_gatk_call.py", args, cwd=OUTPUT_DIR, env=env)

    if STEPS["step4_consensus"]:
        print("\n>>> [Step 4/6] Generating Consensus Sequences...")
        args = [
            "-r",
            str(reference_genome),
            "-v",
            str(vcf_dir),
            "-o",
            str(consensus_dir),
            "--min-gq",
            str(STEP4_MIN_GQ),
            "--min-dp",
            str(STEP4_MIN_DP),
            "-t",
            str(STEP4_PARALLEL_SAMPLES),
        ]
        if not STEP4_USE_IUPAC:
            args.append("--no-iupac")
        run_python_script("4_vcf_to_fasta.py", args, cwd=OUTPUT_DIR, env=env)

    if STEPS["step5_qc"]:
        print("\n>>> [Step 5/6] Quality Control...")
        quality_report = OUTPUT_DIR / "quality_report.txt"
        qc_script = STEP_SCRIPT_DIR / "5_check_non_atcg.py"
        command = f"{sys.executable} {qc_script} | tee {quality_report}"
        run_shell_command(command, cwd=consensus_dir, env=env)

    if STEPS["step6_combine"]:
        print("\n>>> [Step 6/6] Combining Sequences...")
        run_python_script("6_combine_fasta.py", [], cwd=OUTPUT_DIR, env=env)

    duration_minutes = (time.time() - start_time) / 60
    summary_lines = [
        f"reference_genome\t{reference_genome}",
        f"fastq_dir\t{fastq_dir}",
        f"bam_dir\t{bam_dir}",
        f"vcf_dir\t{vcf_dir}",
        f"consensus_dir\t{consensus_dir}",
        f"duration_minutes\t{duration_minutes:.2f}",
    ]
    (OUTPUT_DIR / "pipeline_summary.tsv").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print("\n" + "=" * 60)
    print(f"所有已选步骤运行完成！总耗时: {duration_minutes:.2f} 分钟")
    print("=" * 60)


if __name__ == "__main__":
    main()
