#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INPUT_FASTA_FILE = "result.fasta"
OUTPUT_PREFIX = "bootstrap_replicate"
NUM_REPLICATES = 10


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
    for header, seq in original_sequences[1:]:
        if len(seq) != alignment_length:
            raise ValueError(f"序列长度不一致: {header}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    column_indices = list(range(alignment_length))
    for i in range(1, NUM_REPLICATES + 1):
        bootstrap_indices = random.choices(column_indices, k=alignment_length)
        output_filename = OUTPUT_DIR / f"{OUTPUT_PREFIX}_{i}.fasta"
        with output_filename.open("w", encoding="utf-8") as handle:
            for header, original_seq in original_sequences:
                handle.write(f">{header}\n{''.join(original_seq[j] for j in bootstrap_indices)}\n")
    print("Bootstrap 文件生成完成。")


if __name__ == "__main__":
    main()
