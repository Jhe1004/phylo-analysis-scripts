#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path

from Bio import SeqIO


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
FASTA_FILE = "input.fasta"
MAPPING_FILE = "mapping.csv"
OUTPUT_FASTA_FILE = "output_renamed_linear.fasta"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def create_replacement_map(mapping_file: Path) -> dict[str, str]:
    replacement_dict: dict[str, str] = {}
    with mapping_file.open("r", encoding="utf-8") as infile:
        reader = csv.reader(infile)
        for rows in reader:
            if len(rows) >= 2:
                replacement_dict[rows[0].strip()] = rows[1].strip()
    return replacement_dict


def main() -> None:
    fasta_file = INPUT_DIR / FASTA_FILE
    mapping_file = INPUT_DIR / MAPPING_FILE
    if not fasta_file.exists() or not mapping_file.exists():
        raise FileNotFoundError("输入 FASTA 或映射文件不存在。")
    name_map = create_replacement_map(mapping_file)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / OUTPUT_FASTA_FILE
    records = []
    for record in SeqIO.parse(str(fasta_file), "fasta"):
        if record.id in name_map:
            record.id = name_map[record.id]
            record.description = ""
        records.append(record)
    with output_file.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.id}\n{record.seq}\n")
    print("FASTA 标题重命名完成。")


if __name__ == "__main__":
    main()
