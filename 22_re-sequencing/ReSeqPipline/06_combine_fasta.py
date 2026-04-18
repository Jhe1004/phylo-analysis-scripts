#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import glob
import logging
import os
from collections import OrderedDict
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def ensure_exists(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"未找到{description}: {path}")


# =============================
# 参数配置区
# =============================
CONDA_ENV_NAME = "reseq"
CONSENSUS_DIRECTORY = "output/consensus_fasta_output"


OUTPUT_FILE = "combined_sequences.fasta"


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


LOGGER = setup_logger("06_combine_fasta")


def parse_fasta(filepath: str) -> OrderedDict[str, str] | None:
    sequences: OrderedDict[str, str] = OrderedDict()
    current_id = None
    current_seq_parts: list[str] = []
    try:
        with open(filepath, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_id is not None:
                        sequences[current_id] = "".join(current_seq_parts)
                    current_id = line[1:].split()[0]
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_seq_parts)
    except Exception as exc:
        LOGGER.error(f"读取文件 {os.path.basename(filepath)} 时出错: {exc}")
        return None
    return sequences


def main() -> None:
    consensus_dir = (SCRIPT_DIR / CONSENSUS_DIRECTORY).resolve()
    output_file = (consensus_dir.parent / OUTPUT_FILE).resolve()
    ensure_exists(consensus_dir, "consensus FASTA 目录")

    fasta_files = sorted(glob.glob(os.path.join(str(consensus_dir), "*.fasta")))
    if not fasta_files:
        raise FileNotFoundError(f"在目录 {consensus_dir} 中没有找到 .fasta 文件。")

    all_sample_data: OrderedDict[str, OrderedDict[str, str]] = OrderedDict()
    scaffold_order = None

    for filepath in fasta_files:
        sample_name = os.path.splitext(os.path.basename(filepath))[0]
        if sample_name.endswith("_consensus"):
            sample_name = sample_name[:-len("_consensus")]
        parsed_seqs = parse_fasta(filepath)
        if not parsed_seqs:
            continue
        all_sample_data[sample_name] = parsed_seqs
        if scaffold_order is None:
            scaffold_order = list(parsed_seqs.keys())

    if not all_sample_data or scaffold_order is None:
        raise RuntimeError("没有成功解析任何 FASTA 文件。")

    with open(output_file, "w", encoding="utf-8") as outfile:
        for sample_name, seqs in all_sample_data.items():
            full_sequence = "".join(seqs[scaffold_id] for scaffold_id in scaffold_order if scaffold_id in seqs)
            outfile.write(f">{sample_name}\n")
            for i in range(0, len(full_sequence), 80):
                outfile.write(full_sequence[i:i + 80] + "\n")

    LOGGER.info(f"成功输出合并结果: {output_file}")


if __name__ == "__main__":
    main()
