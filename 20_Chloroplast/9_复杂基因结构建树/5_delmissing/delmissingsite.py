#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from multiprocessing import Pool
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INPUT_EXTENSIONS = [".fasta", ".fas", ".fa"]
MISSING_PROPORTION_THRESHOLD = 0.2
PROCESS_COUNT = 8


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def filter_alignment(input_file: Path) -> tuple[str, int, int]:
    alignment = AlignIO.read(str(input_file), "fasta")
    total_columns = alignment.get_alignment_length()
    kept_indices = []
    for column_index in range(total_columns):
        column = alignment[:, column_index]
        missing_count = sum(base.upper() not in {"A", "T", "C", "G"} for base in column)
        if missing_count / len(column) < MISSING_PROPORTION_THRESHOLD:
            kept_indices.append(column_index)
    output_records = []
    for record in alignment:
        filtered_sequence = "".join(record.seq[index] for index in kept_indices)
        if set(filtered_sequence) <= {"-", "?", "N", "n"}:
            continue
        output_records.append(SeqRecord(Seq(filtered_sequence), id=record.id, description=""))
    output_path = OUTPUT_DIR / f"{input_file.stem}.fas"
    SeqIO.write(output_records, str(output_path), "fasta")
    return input_file.name, total_columns, len(kept_indices)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    input_files = []
    for extension in INPUT_EXTENSIONS:
        input_files.extend(sorted(INPUT_DIR.glob(f"*{extension}")))
    if not input_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到待过滤的比对文件。")

    with Pool(processes=PROCESS_COUNT) as pool:
        results = pool.map(filter_alignment, input_files)

    summary_lines = ["file\toriginal_length\tkept_length"]
    for file_name, original_length, kept_length in results:
        summary_lines.append(f"{file_name}\t{original_length}\t{kept_length}")
    (OUTPUT_DIR / "summary.tsv").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("缺失位点过滤完成。")


if __name__ == "__main__":
    main()
