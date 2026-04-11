#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import concurrent.futures
import glob
import logging
import os
import subprocess

from pipeline_utils import (
    SCRIPT_DIR,
    build_runtime_env,
    ensure_dir,
    ensure_exists,
)


# =============================
# 参数配置区
# =============================
CONDA_ENV_NAME = "reseq"
REFERENCE_GENOME = "input/reference/ref.fasta"
FASTQ_DIRECTORY = "input/fastq"
OUTPUT_DIRECTORY = "output/bam_output"

PLUS_TAG = "_1.clean.fq.gz"
MINUS_TAG = "_2.clean.fq.gz"
THREADS_PER_SAMPLE = 20
PARALLEL_SAMPLES = 10


def setup_logger(script_name: str) -> logging.Logger:
    log_file = SCRIPT_DIR / f"{script_name}.log"
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    logger.info(f"日志文件: {log_file}")
    return logger


LOGGER = setup_logger("02_bwa_mapping")
RUNTIME_ENV = build_runtime_env(CONDA_ENV_NAME)


def run_command(command: str, description: str, sample_name: str) -> None:
    prefix = f"[{sample_name}]"
    LOGGER.info(f"\n{prefix} [子步骤] {description}")
    LOGGER.debug(f"{prefix} [命令] {command}")
    result = subprocess.run(command, shell=True, cwd=str(SCRIPT_DIR), env=RUNTIME_ENV)
    if result.returncode != 0:
        raise RuntimeError(f"{prefix} 失败: {description}")


def find_samples(folder: str, plus_tag: str, minus_tag: str) -> list[dict[str, str]]:
    samples_dict: dict[str, dict[str, str]] = {}
    tag_pairs = [
        (plus_tag, minus_tag),
        ("_1.clean.fq.gz", "_2.clean.fq.gz"),
        ("_1.fq.gz", "_2.fq.gz"),
        ("_R1.fastq.gz", "_R2.fastq.gz"),
        ("_R1.fq.gz", "_R2.fq.gz"),
    ]

    seen_pairs = set()
    unique_tag_pairs = []
    for pair in tag_pairs:
        if pair not in seen_pairs:
            unique_tag_pairs.append(pair)
            seen_pairs.add(pair)

    for p_tag, m_tag in unique_tag_pairs:
        for filepath in glob.glob(os.path.join(folder, f"*{p_tag}")):
            filename = os.path.basename(filepath)
            sample_name = filename.replace(p_tag, "")
            if sample_name in samples_dict:
                continue
            minus_file = os.path.join(folder, f"{sample_name}{m_tag}")
            if os.path.exists(minus_file):
                samples_dict[sample_name] = {
                    "name": sample_name,
                    "r1": filepath,
                    "r2": minus_file,
                }

    LOGGER.info(f"匹配到 {len(samples_dict)} 个样本。")
    return list(samples_dict.values())


def process_sample(sample: dict[str, str], reference: str, threads: int, output_folder: str) -> str | None:
    sample_name = sample["name"]
    raw_bam = os.path.join(output_folder, f"{sample_name}.raw.bam")
    sorted_bam = os.path.join(output_folder, f"{sample_name}.sorted.bam")
    final_bam = os.path.join(output_folder, f"{sample_name}.bam")

    if os.path.exists(final_bam) and os.path.getsize(final_bam) > 0:
        LOGGER.info(f"[{sample_name}] 最终 BAM 已存在，跳过。")
        return final_bam

    try:
        rg_string = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tLB:WGS\\tPL:ILLUMINA"
        bwa_cmd = (
            f"bwa mem -t {threads} -M -R '{rg_string}' "
            f"{reference} {sample['r1']} {sample['r2']} | "
            f"samtools view -bS - > {raw_bam}"
        )
        run_command(bwa_cmd, "BWA MEM 比对", sample_name)
        run_command(f"samtools sort -@ {threads} -o {sorted_bam} {raw_bam}", "Samtools 排序", sample_name)

        if os.path.exists(raw_bam):
            os.remove(raw_bam)
        if os.path.exists(sorted_bam):
            os.rename(sorted_bam, final_bam)

        run_command(f"samtools index {final_bam}", "Samtools 索引", sample_name)
        LOGGER.info(f"[{sample_name}] [完成] 输出: {final_bam}")
        return final_bam
    except Exception as exc:
        LOGGER.error(f"[{sample_name}] [严重错误] 样本处理中断: {exc}")
        return None


def main() -> None:
    reference_genome = (SCRIPT_DIR / REFERENCE_GENOME).resolve()
    fastq_dir = (SCRIPT_DIR / FASTQ_DIRECTORY).resolve()
    output_dir = (SCRIPT_DIR / OUTPUT_DIRECTORY).resolve()

    ensure_exists(reference_genome, "参考基因组文件")
    ensure_exists(fastq_dir, "FASTQ 输入目录")
    ensure_dir(output_dir)

    samples = find_samples(str(fastq_dir), PLUS_TAG, MINUS_TAG)
    if not samples:
        raise FileNotFoundError(f"未找到匹配的 FASTQ 样本: {fastq_dir}")

    num_parallel = min(PARALLEL_SAMPLES, len(samples))
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel) as executor:
        futures = {
            executor.submit(process_sample, sample, str(reference_genome), THREADS_PER_SAMPLE, str(output_dir)): sample["name"]
            for sample in samples
        }
        for future in concurrent.futures.as_completed(futures):
            future.result()


if __name__ == "__main__":
    main()
