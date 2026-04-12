#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import glob
import logging
import os
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEPENDENCIES_DIR = SCRIPT_DIR / "dependencies"


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


def build_runtime_env(conda_env_name: str) -> dict[str, str]:
    env = os.environ.copy()
    path_parts = []

    local_bin = DEPENDENCIES_DIR / "bin"
    if local_bin.exists():
        path_parts.append(str(local_bin))

    conda_bin = resolve_conda_bin(conda_env_name)
    if conda_bin is not None:
        path_parts.append(str(conda_bin))

    path_parts.append(env.get("PATH", ""))
    env["PATH"] = ":".join(part for part in path_parts if part)
    return env


def ensure_exists(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"未找到{description}: {path}")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def ensure_tools_available(tools: list[str]) -> None:
    missing_tools = []
    runtime_path = RUNTIME_ENV.get("PATH", os.environ.get("PATH", ""))
    for tool in tools:
        if shutil.which(tool, path=runtime_path) is None:
            missing_tools.append(tool)
    if missing_tools:
        raise EnvironmentError(
            "缺少以下依赖软件，请先在 conda 环境中安装: "
            + ", ".join(missing_tools)
        )


# =============================
# 参数配置区
# =============================
CONDA_ENV_NAME = "reseq"
REFERENCE_GENOME = "ref.fasta"
BAM_DIRECTORY = "bam_output"
OUTPUT_DIRECTORY = "vcf_output"

PARALLEL_SAMPLES = 20
GATK_THREADS = 6
JOINT_CALLING = False
SINGLE_VCF = False
COMPRESS_VCF_OUTPUT = False


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


def format_progress(done: int, total: int) -> str:
    if total <= 0:
        return "0/0 (0.0%)"
    return f"{done}/{total} ({done / total:.1%})"


def run_command(command: str, description: str) -> None:
    result = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        env=RUNTIME_ENV,
    )
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", errors="replace")
        raise RuntimeError(
            f"{description} 失败\n"
            f"[具体命令] {command}\n"
            f"[Stderr] {stderr}"
        )


def get_chromosomes(fai_path: str) -> list[str]:
    chroms = []
    with open(fai_path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if parts:
                chroms.append(parts[0])
    return chroms


def get_max_contig_length(fai_path: str) -> int:
    max_length = 0
    with open(fai_path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                max_length = max(max_length, int(parts[1]))
    return max_length


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
    suffix = ".g.vcf.gz" if emit_gvcf and COMPRESS_VCF_OUTPUT else ".g.vcf" if emit_gvcf else ".vcf.gz" if COMPRESS_VCF_OUTPUT else ".vcf"
    output_file = os.path.join(shard_dir, f"{sample_name}.{chrom}{suffix}")
    if os.path.exists(output_file) and os.path.getsize(output_file) > 1000:
        return output_file

    erc_option = "--emit-ref-confidence GVCF" if emit_gvcf else ""
    cmd = (
        f"gatk HaplotypeCaller -R {shlex.quote(reference)} "
        f"-I {shlex.quote(bam_path)} -O {shlex.quote(output_file)} "
        f"-L {shlex.quote(chrom)} {erc_option} --native-pair-hmm-threads {gatk_threads}"
    )
    run_command(cmd, f"{sample_name} - {chrom}")
    return output_file


def gather_shards(sample_name: str, shards: list[str], output_folder: str, emit_gvcf: bool) -> str:
    suffix = ".g.vcf.gz" if emit_gvcf and COMPRESS_VCF_OUTPUT else ".g.vcf" if emit_gvcf else ".vcf.gz" if COMPRESS_VCF_OUTPUT else ".vcf"
    final_output = os.path.join(output_folder, f"{sample_name}{suffix}")
    inputs = " ".join([f"-I {shlex.quote(shard)}" for shard in shards])
    cmd = f"gatk GatherVcfsCloud {inputs} -O {shlex.quote(final_output)} --ignore-safety-checks"
    if subprocess.run(
        "gatk GatherVcfsCloud --help",
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env=RUNTIME_ENV,
    ).returncode != 0:
        cmd = cmd.replace("GatherVcfsCloud", "GatherVcfs")
    run_command(cmd, f"合并样本 {sample_name}")
    run_command(f"gatk IndexFeatureFile -I {shlex.quote(final_output)}", f"建立索引 {sample_name}")
    return final_output


def combine_gvcfs(gvcf_files: list[str], reference: str, output_folder: str) -> str:
    combined_gvcf = os.path.join(output_folder, "cohort.g.vcf.gz" if COMPRESS_VCF_OUTPUT else "cohort.g.vcf")
    list_file = os.path.join(output_folder, "gvcfs.list")
    with open(list_file, "w", encoding="utf-8") as handle:
        for gvcf in gvcf_files:
            handle.write(gvcf + "\n")
    run_command(
        f"gatk CombineGVCFs -R {shlex.quote(reference)} "
        f"-V {shlex.quote(list_file)} -O {shlex.quote(combined_gvcf)}",
        "GATK CombineGVCFs",
    )
    return combined_gvcf


def genotype_gvcfs(combined_gvcf: str, reference: str, output_folder: str) -> str:
    final_vcf = os.path.join(output_folder, "all_samples.vcf.gz" if COMPRESS_VCF_OUTPUT else "all_samples.vcf")
    run_command(
        f"gatk GenotypeGVCFs -R {shlex.quote(reference)} "
        f"-V {shlex.quote(combined_gvcf)} -O {shlex.quote(final_vcf)} "
        f"--include-non-variant-sites",
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
    ensure_tools_available(["gatk", "java", "samtools", "bgzip", "tabix"])

    fai_file = f"{reference_genome}.fai"
    dict_file = reference_genome.with_suffix(".dict")
    ensure_exists(Path(fai_file), ".fai 索引文件")
    ensure_exists(dict_file, ".dict 文件")

    chroms = get_chromosomes(fai_file)
    max_contig_length = get_max_contig_length(fai_file)
    bam_files = sorted(glob.glob(os.path.join(str(bam_dir), "*.bam")))
    bam_files = [bam for bam in bam_files if not any(tag in bam for tag in [".raw.", ".sorted."])]
    if not bam_files:
        raise FileNotFoundError(f"未找到 BAM 文件: {bam_dir}")

    LOGGER.info(
        "运行参数: BAM=%s, 染色体=%s, 最长染色体=%s, PARALLEL_SAMPLES=%s, GATK_THREADS=%s, JOINT_CALLING=%s, SINGLE_VCF=%s, COMPRESS_VCF_OUTPUT=%s",
        len(bam_files),
        len(chroms),
        max_contig_length,
        PARALLEL_SAMPLES,
        GATK_THREADS,
        JOINT_CALLING,
        SINGLE_VCF,
        COMPRESS_VCF_OUTPUT,
    )

    if COMPRESS_VCF_OUTPUT and max_contig_length > 512_000_000:
        raise ValueError(
            f"参考基因组最长染色体为 {max_contig_length} bp，超过 tabix/TBI 对压缩 VCF 常见支持范围。"
            "请将 COMPRESS_VCF_OUTPUT 设为 False，输出未压缩 VCF/GVCF 及 .idx 索引。"
        )

    emit_gvcf = JOINT_CALLING or (not SINGLE_VCF)
    all_tasks = [(bam, chrom) for bam in bam_files for chrom in chroms]
    shard_results: dict[str, dict[str, str]] = {}
    sample_names = [os.path.basename(bam).replace(".bam", "") for bam in bam_files]

    LOGGER.info("检测到样本: %s", ", ".join(sample_names))
    LOGGER.info(
        "分片任务总数: %s，预计每个样本 %s 个染色体分片",
        len(all_tasks),
        len(chroms),
    )
    LOGGER.info("开始并行运行 HaplotypeCaller 分片任务")

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
        completed_tasks = 0
        progress_interval = max(1, min(25, len(all_tasks) // 20 or 1))
        for future in as_completed(futures):
            bam, chrom = futures[future]
            out_file = future.result()
            shard_results.setdefault(bam, {})[chrom] = out_file
            completed_tasks += 1

            sample_name = os.path.basename(bam).replace(".bam", "")
            sample_done = len(shard_results[bam])
            if completed_tasks == 1 or completed_tasks % progress_interval == 0 or completed_tasks == len(all_tasks):
                LOGGER.info(
                    "分片进度: %s；刚完成 %s:%s；该样本进度 %s",
                    format_progress(completed_tasks, len(all_tasks)),
                    sample_name,
                    chrom,
                    format_progress(sample_done, len(chroms)),
                )

    final_sample_files = []
    LOGGER.info("所有分片任务完成，开始按样本合并分片")
    for index, bam in enumerate(bam_files, start=1):
        sample_name = os.path.basename(bam).replace(".bam", "")
        sample_shards = [shard_results[bam][chrom] for chrom in chroms]
        LOGGER.info(
            "开始合并样本 %s，样本进度 %s",
            sample_name,
            format_progress(index, len(bam_files)),
        )
        final_sample_files.append(gather_shards(sample_name, sample_shards, str(output_dir), emit_gvcf))
        LOGGER.info(
            "样本合并完成: %s，输出文件: %s",
            sample_name,
            final_sample_files[-1],
        )

    if JOINT_CALLING and len(final_sample_files) > 1:
        LOGGER.info("开始联合分型：CombineGVCFs")
        combined_gvcf = combine_gvcfs(final_sample_files, str(reference_genome), str(output_dir))
        LOGGER.info("CombineGVCFs 完成: %s", combined_gvcf)
        LOGGER.info("开始联合分型：GenotypeGVCFs")
        genotype_gvcfs(combined_gvcf, str(reference_genome), str(output_dir))
        LOGGER.info("GenotypeGVCFs 完成")

    LOGGER.info(
        "流程完成: 成功生成 %s 个样本结果文件，输出目录: %s",
        len(final_sample_files),
        output_dir,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception:
        LOGGER.exception("流程异常终止")
        raise
