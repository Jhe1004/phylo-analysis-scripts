#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from collections import defaultdict
from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
FASTA_EXTENSION = ".fasta"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def extract_sequences_to_files(folder_path: Path, output_folder: Path) -> None:
    sequences: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for file in sorted(folder_path.glob(f"*{FASTA_EXTENSION}")):
        with file.open("r", encoding="utf-8") as fasta:
            seq_name = ""
            for line in fasta:
                if line.startswith(">"):
                    seq_name = line.strip().replace(">", "")
                else:
                    sequences[seq_name].append((file.name, line.strip()))

    output_folder.mkdir(parents=True, exist_ok=True)
    for seq_name, seq_list in sequences.items():
        with (output_folder / f"{seq_name}.fasta").open("w", encoding="utf-8") as output:
            for source_file, sequence in seq_list:
                output.write(f">{source_file}\n{sequence}\n")


def main() -> None:
    fasta_files = sorted(INPUT_DIR.glob(f"*{FASTA_EXTENSION}"))
    if not fasta_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 *{FASTA_EXTENSION} 文件。")
    extract_sequences_to_files(INPUT_DIR, OUTPUT_DIR)
    print("已按基因名汇总序列。")


if __name__ == "__main__":
    main()
