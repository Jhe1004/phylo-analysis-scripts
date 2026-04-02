#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from Bio import SeqIO


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INPUT_EXTENSIONS = [".fasta", ".fa", ".cds"]


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def calculate_n50(lengths: list[int]) -> int | None:
    lengths_sorted = sorted(lengths, reverse=True)
    half_total_length = sum(lengths) / 2
    running_sum = 0
    for length in lengths_sorted:
        running_sum += length
        if running_sum >= half_total_length:
            return length
    return None


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = ["file\tsequence_count\tn50"]
    fasta_files = []
    for ext in INPUT_EXTENSIONS:
        fasta_files.extend(sorted(INPUT_DIR.glob(f"*{ext}")))
    if not fasta_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 FASTA 文件。")
    for fasta_path in fasta_files:
        seq_lengths = [len(record.seq) for record in SeqIO.parse(str(fasta_path), "fasta")]
        results.append(f"{fasta_path.name}\t{len(seq_lengths)}\t{calculate_n50(seq_lengths)}")
    (OUTPUT_DIR / "assembly_stats.tsv").write_text("\n".join(results) + "\n", encoding="utf-8")
    print("N50 统计完成。")


if __name__ == "__main__":
    main()
