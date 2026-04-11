#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import logging
import os
import subprocess
import sys

from pipeline_utils import SCRIPT_DIR, build_runtime_env, ensure_exists


# =============================
# 参数配置区
# 直接修改这里即可，无需命令行传参
# =============================
CONDA_ENV_NAME = "reseq"
REFERENCE_GENOME = "input/reference/ref.fasta"
BWA_INDEX_ALGORITHM = "bwtsw"


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


LOGGER = setup_logger("01_index_reference")


def run_command(command: str, description: str, env: dict[str, str]) -> None:
    LOGGER.info(f"\n{'=' * 60}")
    LOGGER.info(f"[步骤] {description}")
    LOGGER.debug(f"[命令] {command}")
    LOGGER.info("=" * 60)
    result = subprocess.run(command, shell=True, cwd=str(SCRIPT_DIR), env=env)
    if result.returncode != 0:
        LOGGER.error(f"[错误] 命令执行失败，返回码: {result.returncode}")
        sys.exit(1)
    LOGGER.info(f"[完成] {description}")


def main() -> None:
    reference_genome = (SCRIPT_DIR / REFERENCE_GENOME).resolve()
    ensure_exists(reference_genome, "参考基因组文件")
    env = build_runtime_env(CONDA_ENV_NAME)

    ref_dir = reference_genome.parent
    ref_name_no_ext = os.path.splitext(reference_genome.name)[0]
    dict_path = ref_dir / f"{ref_name_no_ext}.dict"

    LOGGER.info(f"\n参考基因组: {reference_genome}")
    LOGGER.info(f"索引文件将保存在: {ref_dir}")

    run_command(
        f"bwa index -a {BWA_INDEX_ALGORITHM} {reference_genome}",
        "创建 BWA 索引",
        env,
    )
    run_command(
        f"samtools faidx {reference_genome}",
        "创建 Samtools FASTA 索引 (.fai)",
        env,
    )

    if dict_path.exists():
        LOGGER.info(f"[信息] 发现已存在的 .dict 文件，将删除后重新生成: {dict_path}")
        dict_path.unlink()

    gatk_cmd = f"gatk CreateSequenceDictionary -R {reference_genome} -O {dict_path}"
    result = subprocess.run(gatk_cmd, shell=True, cwd=str(SCRIPT_DIR), env=env)
    if result.returncode != 0:
        LOGGER.info("[信息] GATK 命令失败，尝试使用 Picard...")
        run_command(
            f"picard CreateSequenceDictionary R={reference_genome} O={dict_path}",
            "创建 Picard 序列字典 (.dict)",
            env,
        )
    else:
        LOGGER.info("[完成] 创建 GATK 序列字典 (.dict)")

    LOGGER.info("\n" + "=" * 60)
    LOGGER.info("所有索引文件创建完成！")
    LOGGER.info("=" * 60)


if __name__ == "__main__":
    main()
