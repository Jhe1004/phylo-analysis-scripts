#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import logging
import os
from collections import Counter
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def ensure_exists(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"未找到{description}: {path}")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


# =============================
# 参数配置区
# =============================
CONDA_ENV_NAME = "reseq"
CONSENSUS_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/22_re-sequencing/ReSeqPipline/output/consensus_fasta_output"
REPORT_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/22_re-sequencing/ReSeqPipline/output/quality_report.txt"


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


LOGGER = setup_logger("05_check_non_atcg")


def is_fasta(filename: str) -> bool:
    return filename.lower().endswith((".fasta", ".fa", ".fna", ".fas"))


def count_non_atcg(filepath: str) -> tuple[int, int] | None:
    total_bases = 0
    base_counts: Counter[str] = Counter()
    try:
        with open(filepath, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if line.startswith(">") or not line:
                    continue
                seq_upper = line.upper()
                total_bases += len(seq_upper)
                base_counts.update(seq_upper)
    except Exception as exc:
        LOGGER.error(f"读取文件 {os.path.basename(filepath)} 时出错: {exc}")
        return None

    atcg_count = base_counts["A"] + base_counts["T"] + base_counts["C"] + base_counts["G"]
    non_atcg_count = total_bases - atcg_count
    return total_bases, non_atcg_count


def main() -> None:
    consensus_dir = (SCRIPT_DIR / CONSENSUS_DIRECTORY).resolve()
    report_file = (SCRIPT_DIR / REPORT_FILE).resolve()

    ensure_exists(consensus_dir, "consensus FASTA 目录")
    ensure_dir(report_file.parent)

    fasta_files = [name for name in os.listdir(consensus_dir) if (consensus_dir / name).is_file() and is_fasta(name)]
    fasta_files.sort()
    if not fasta_files:
        LOGGER.info("当前文件夹中没有找到 FASTA 文件。")
        return

    lines = []
    header = ["Filename", "Total Bases", "Non-ATCG", "Ratio (%)"]
    table_data = []
    grand_total_bases = 0
    grand_total_non_atcg = 0

    for file in fasta_files:
        result = count_non_atcg(str(consensus_dir / file))
        if result is None:
            continue
        total, non_atcg = result
        grand_total_bases += total
        grand_total_non_atcg += non_atcg
        ratio = (non_atcg / total * 100) if total > 0 else 0
        table_data.append([file, total, non_atcg, ratio])

    col_widths = [len(h) for h in header]
    for row in table_data:
        col_widths[0] = max(col_widths[0], len(str(row[0])))
        col_widths[1] = max(col_widths[1], len(f"{row[1]:,}"))
        col_widths[2] = max(col_widths[2], len(f"{row[2]:,}"))
        col_widths[3] = max(col_widths[3], len(f"{row[3]:.4f}"))

    row_format = " | ".join(
        [f"{{:<{col_widths[0]}}}", f"{{:>{col_widths[1]}}}", f"{{:>{col_widths[2]}}}", f"{{:>{col_widths[3]}}}"]
    )
    separator = "-" * (sum(col_widths) + 9)
    lines.append(separator)
    lines.append(row_format.format(*header))
    lines.append(separator)
    for row in table_data:
        lines.append(row_format.format(row[0], f"{row[1]:,}", f"{row[2]:,}", f"{row[3]:.4f}"))
    grand_ratio = (grand_total_non_atcg / grand_total_bases * 100) if grand_total_bases > 0 else 0
    lines.append(separator)
    lines.append(row_format.format("TOTAL SUMMARY", f"{grand_total_bases:,}", f"{grand_total_non_atcg:,}", f"{grand_ratio:.4f}"))
    lines.append(separator)

    report_text = "\n".join(lines) + "\n"
    report_file.write_text(report_text, encoding="utf-8")
    print(report_text, end="")
    LOGGER.info(f"统计结果已写入: {report_file}")


if __name__ == "__main__":
    main()
