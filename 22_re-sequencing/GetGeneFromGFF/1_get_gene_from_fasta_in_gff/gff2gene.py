#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INPUT_GFF_FILE = "ref.gff"
FASTA_EXTENSION = ".fasta"
OUTPUT_PREFIX = "extracted_"
MAX_GENE_LENGTH = 10000
TARGET_FEATURE = "mRNA"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def extract_gene_positions(gff_file: Path) -> list[tuple[int, int]]:
    gene_positions = []
    with gff_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[2] != TARGET_FEATURE:
                continue
            start, end = int(fields[3]), int(fields[4])
            if end - start + 1 <= MAX_GENE_LENGTH:
                gene_positions.append((start, end))
    return gene_positions


def extract_genes_from_fasta(gene_positions: list[tuple[int, int]], fasta_file: Path, output_file: Path) -> None:
    with fasta_file.open("r", encoding="utf-8") as fasta, output_file.open("w", encoding="utf-8") as output:
        fasta.readline()
        sequence = fasta.read().replace("\n", "")
        for index, (start, end) in enumerate(gene_positions, start=1):
            gene_sequence = sequence[start - 1 : end]
            output.write(f">gene_{index}\n{gene_sequence}\n")


def main() -> None:
    gff_file = INPUT_DIR / INPUT_GFF_FILE
    if not gff_file.exists():
        raise FileNotFoundError(f"未找到输入 GFF 文件: {gff_file}")

    fasta_files = sorted(INPUT_DIR.glob(f"*{FASTA_EXTENSION}"))
    if not fasta_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 *{FASTA_EXTENSION} 文件。")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    gene_positions = extract_gene_positions(gff_file)
    for fasta_file in fasta_files:
        output_file = OUTPUT_DIR / f"{OUTPUT_PREFIX}{fasta_file.name}"
        extract_genes_from_fasta(gene_positions, fasta_file, output_file)
        print(f"已处理: {fasta_file.name}")

    print("所有 FASTA 文件提取完成。")


if __name__ == "__main__":
    main()
