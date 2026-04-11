#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import glob
import logging
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

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
BAM_DIRECTORY = "output/bam_output"
OUTPUT_DIRECTORY = "output/vcf_output"

PARALLEL_SAMPLES = 20
GATK_THREADS = 6
JOINT_CALLING = False
SINGLE_VCF = False


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


LOGGER = setup_logger("03_gatk_call")
RUNTIME_ENV = build_runtime_env(CONDA_ENV_NAME)


def run_command(command: str, description: str) -> None:
    result = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        env=RUNTIME_ENV,
    )
    if result.returncode != 0:
        LOGGER.error(f"\n[错误] {description} 失败")
        LOGGER.error(f"[具体命令] {command}")
        LOGGER.error(f"[Stderr] {result.stderr.decode('utf-8')}")
        sys.exit(1)


def get_chromosomes(fai_path: str) -> list[str]:
    chroms = []
    with open(fai_path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if parts:
                chroms.append(parts[0])
    return chroms


def call_variant_shard(
    bam_path: str,
    reference: str,
    output_folder: str,
    chrom: str,
    emit_gvcf: bool,
    gatk_threads: int,
) -> str:
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    shard_dir = os.path.join(output_folder, "shards", sample_name)
    os.makedirs(shard_dir, exist_ok=True)
    suffix = ".g.vcf.gz" if emit_gvcf else ".vcf.gz"
    output_file = os.path.join(shard_dir, f"{sample_name}.{chrom}{suffix}")
    if os.path.exists(output_file) and os.path.getsize(output_file) > 1000:
        return output_file

    erc_option = "--emit-ref-confidence GVCF" if emit_gvcf else ""
    cmd = (
        f"gatk HaplotypeCaller -R {reference} -I {bam_path} -O {output_file} "
        f"-L {chrom} {erc_option} --native-pair-hmm-threads {gatk_threads}"
    )
    run_command(cmd, f"{sample_name} - {chrom}")
    return output_file


def gather_shards(sample_name: str, shards: list[str], output_folder: str, emit_gvcf: bool) -> str:
    suffix = ".g.vcf.gz" if emit_gvcf else ".vcf.gz"
    final_output = os.path.join(output_folder, f"{sample_name}{suffix}")
    inputs = " ".join([f"-I {shard}" for shard in shards])
    cmd = f"gatk GatherVcfsCloud {inputs} -O {final_output} --ignore-safety-checks"
    if subprocess.run(
        "gatk GatherVcfsCloud --help",
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env=RUNTIME_ENV,
    ).returncode != 0:
        cmd = cmd.replace("GatherVcfsCloud", "GatherVcfs")
    run_command(cmd, f"合并样本 {sample_name}")
    run_command(f"gatk IndexFeatureFile -I {final_output}", f"建立索引 {sample_name}")
    return final_output


def combine_gvcfs(gvcf_files: list[str], reference: str, output_folder: str) -> str:
    combined_gvcf = os.path.join(output_folder, "cohort.g.vcf.gz")
    list_file = os.path.join(output_folder, "gvcfs.list")
    with open(list_file, "w", encoding="utf-8") as handle:
        for gvcf in gvcf_files:
            handle.write(gvcf + "\n")
    run_command(
        f"gatk CombineGVCFs -R {reference} -V {list_file} -O {combined_gvcf}",
        "GATK CombineGVCFs",
    )
    return combined_gvcf


def genotype_gvcfs(combined_gvcf: str, reference: str, output_folder: str) -> str:
    final_vcf = os.path.join(output_folder, "all_samples.vcf.gz")
    run_command(
        f"gatk GenotypeGVCFs -R {reference} -V {combined_gvcf} -O {final_vcf} --include-non-variant-sites",
        "GATK GenotypeGVCFs",
    )
    return final_vcf


def main() -> None:
    reference_genome = (SCRIPT_DIR / REFERENCE_GENOME).resolve()
    bam_dir = (SCRIPT_DIR / BAM_DIRECTORY).resolve()
    output_dir = (SCRIPT_DIR / OUTPUT_DIRECTORY).resolve()

    ensure_exists(reference_genome, "参考基因组文件")
    ensure_exists(bam_dir, "BAM 输入目录")
    ensure_dir(output_dir)

    fai_file = f"{reference_genome}.fai"
    dict_file = reference_genome.with_suffix(".dict")
    ensure_exists(Path(fai_file), ".fai 索引文件")
    ensure_exists(dict_file, ".dict 文件")

    chroms = get_chromosomes(fai_file)
    bam_files = sorted(glob.glob(os.path.join(str(bam_dir), "*.bam")))
    bam_files = [bam for bam in bam_files if not any(tag in bam for tag in [".raw.", ".sorted."])]
    if not bam_files:
        raise FileNotFoundError(f"未找到 BAM 文件: {bam_dir}")

    emit_gvcf = JOINT_CALLING or (not SINGLE_VCF)
    all_tasks = [(bam, chrom) for bam in bam_files for chrom in chroms]
    shard_results: dict[str, dict[str, str]] = {}

    with ProcessPoolExecutor(max_workers=PARALLEL_SAMPLES) as executor:
        futures = {
            executor.submit(
                call_variant_shard,
                bam,
                str(reference_genome),
                str(output_dir),
                chrom,
                emit_gvcf,
                GATK_THREADS,
            ): (bam, chrom)
            for bam, chrom in all_tasks
        }
        for future in as_completed(futures):
            bam, chrom = futures[future]
            out_file = future.result()
            shard_results.setdefault(bam, {})[chrom] = out_file

    final_sample_files = []
    for bam in bam_files:
        sample_name = os.path.basename(bam).replace(".bam", "")
        sample_shards = [shard_results[bam][chrom] for chrom in chroms]
        final_sample_files.append(gather_shards(sample_name, sample_shards, str(output_dir), emit_gvcf))

    if JOINT_CALLING and len(final_sample_files) > 1:
        combined_gvcf = combine_gvcfs(final_sample_files, str(reference_genome), str(output_dir))
        genotype_gvcfs(combined_gvcf, str(reference_genome), str(output_dir))


if __name__ == "__main__":
    main()
