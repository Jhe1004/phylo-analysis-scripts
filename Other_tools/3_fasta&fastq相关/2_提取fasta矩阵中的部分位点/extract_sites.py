#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
from pathlib import Path


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/3_fasta&fastq相关/2_提取fasta矩阵中的部分位点/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/3_fasta&fastq相关/2_提取fasta矩阵中的部分位点/output"
INPUT_FASTA_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/3_fasta&fastq相关/2_提取fasta矩阵中的部分位点/input/result.fasta"
OUTPUT_FASTA_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/3_fasta&fastq相关/2_提取fasta矩阵中的部分位点/output/extracted_sites.fasta"
PERCENTAGE_TO_EXTRACT = 10.0


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def parse_fasta(filename: Path) -> list[tuple[str, str]]:
    sequences = []
    current_seq_parts = []
    current_header = None
    with filename.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_seq_parts)))
                current_header = line[1:]
                current_seq_parts = []
            else:
                current_seq_parts.append(line)
    if current_header is not None:
        sequences.append((current_header, "".join(current_seq_parts)))
    return sequences


def main() -> None:
    input_file = INPUT_DIR / INPUT_FASTA_FILE
    if not input_file.exists():
        raise FileNotFoundError(f"未找到输入 FASTA 文件: {input_file}")
    original_sequences = parse_fasta(input_file)
    alignment_length = len(original_sequences[0][1])
    num_sites_to_extract = int(alignment_length * (PERCENTAGE_TO_EXTRACT / 100.0))
    extracted_indices = sorted(random.sample(range(alignment_length), k=num_sites_to_extract))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / OUTPUT_FASTA_FILE
    with output_file.open("w", encoding="utf-8") as handle:
        for header, original_seq in original_sequences:
            handle.write(f">{header}\n{''.join(original_seq[i] for i in extracted_indices)}\n")
    print("位点提取完成。")


if __name__ == "__main__":
    main()
