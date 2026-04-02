#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
FASTA_EXTENSION = ".fasta"
MAX_MISSING_RATIO = 0.8
MIN_SEQUENCE_COUNT = 4


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def filter_one_fasta(file_path: Path, output_path: Path) -> bool:
    filtered_sequences = []
    with file_path.open("r", encoding="utf-8") as fasta:
        seq_name = ""
        for line in fasta:
            if line.startswith(">"):
                seq_name = line.strip()
            else:
                seq = line.strip()
                if seq and seq.count("?") / len(seq) <= MAX_MISSING_RATIO:
                    filtered_sequences.append((seq_name, seq))

    if len(filtered_sequences) < MIN_SEQUENCE_COUNT:
        return False

    with output_path.open("w", encoding="utf-8") as output:
        for seq_name, seq in filtered_sequences:
            output.write(f"{seq_name}\n{seq}\n")
    return True


def main() -> None:
    fasta_files = sorted(INPUT_DIR.glob(f"*{FASTA_EXTENSION}"))
    if not fasta_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 *{FASTA_EXTENSION} 文件。")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    kept_files = 0
    removed_files = 0
    summary_lines = ["file\tstatus\tkept_sequence_count"]

    for fasta_file in fasta_files:
        output_file = OUTPUT_DIR / fasta_file.name
        kept = filter_one_fasta(fasta_file, output_file)
        if kept:
            kept_files += 1
            kept_count = sum(1 for line in output_file.open("r", encoding="utf-8") if line.startswith(">"))
            summary_lines.append(f"{fasta_file.name}\tkept\t{kept_count}")
        else:
            removed_files += 1
            if output_file.exists():
                output_file.unlink()
            summary_lines.append(f"{fasta_file.name}\tremoved\t0")

    (OUTPUT_DIR / "summary.tsv").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(f"过滤完成，保留 {kept_files} 个文件，移除 {removed_files} 个文件。")


if __name__ == "__main__":
    main()
