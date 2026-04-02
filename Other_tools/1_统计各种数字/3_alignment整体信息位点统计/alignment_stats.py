#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import Counter
from pathlib import Path

from Bio import AlignIO


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INPUT_FASTA_FILE = "alignment.fasta"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def calculate_alignment_stats(file_path: Path) -> dict[str, int]:
    alignment = AlignIO.read(str(file_path), "fasta")
    alignment_length = alignment.get_alignment_length()
    variable_sites = 0
    parsimony_informative_sites = 0
    for i in range(alignment_length):
        column = "".join(base for base in alignment[:, i].upper() if base in "ATCG")
        if not column:
            continue
        counts = Counter(column)
        unique_bases = [base for base in counts if base in "ATCG"]
        if len(unique_bases) > 1:
            variable_sites += 1
            if sum(1 for base in unique_bases if counts[base] >= 2) >= 2:
                parsimony_informative_sites += 1
    return {
        "characters": alignment_length,
        "variable_sites": variable_sites,
        "parsimony_informative_sites": parsimony_informative_sites,
    }


def main() -> None:
    input_file = INPUT_DIR / INPUT_FASTA_FILE
    if not input_file.exists():
        raise FileNotFoundError(f"未找到比对文件: {input_file}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stats = calculate_alignment_stats(input_file)
    lines = [f"{key}\t{value}" for key, value in stats.items()]
    (OUTPUT_DIR / "alignment_stats.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("比对统计完成。")


if __name__ == "__main__":
    main()
