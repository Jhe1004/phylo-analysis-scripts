#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INPUT_FASTA_FILE = "combined_sequences.fas"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    file_path = INPUT_DIR / INPUT_FASTA_FILE
    if not file_path.exists():
        raise FileNotFoundError(f"未找到 FASTA 文件: {file_path}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_output_path = OUTPUT_DIR / f"{file_path.stem}_analysis.csv"
    standard_bases = {"A", "T", "C", "G", "a", "t", "c", "g"}
    sequences = {}
    current_header = None
    with file_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_header = line[1:]
                sequences[current_header] = ""
            elif current_header:
                sequences[current_header] += line
    with csv_output_path.open("w", newline="", encoding="utf-8-sig") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["序列ID", "总长度(bp)", "非'ATCG'字符数", "所占比例(%)"])
        for header, sequence in sequences.items():
            total_length = len(sequence)
            non_standard_count = sum(1 for base in sequence if base not in standard_bases)
            ratio = (non_standard_count / total_length * 100) if total_length else 0
            writer.writerow([header, total_length, non_standard_count, f"{ratio:.2f}"])
    print("gap 信息统计完成。")


if __name__ == "__main__":
    main()
