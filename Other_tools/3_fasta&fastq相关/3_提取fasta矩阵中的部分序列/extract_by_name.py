#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/3_fasta&fastq相关/3_提取fasta矩阵中的部分序列/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/3_fasta&fastq相关/3_提取fasta矩阵中的部分序列/output"
FASTA_INPUT_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/3_fasta&fastq相关/3_提取fasta矩阵中的部分序列/input/all_sequences.fasta"
NAMES_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/3_fasta&fastq相关/3_提取fasta矩阵中的部分序列/input/names_to_extract.txt"


FASTA_OUTPUT_FILE = "extracted_sequences.fasta"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def read_names_to_extract(filename: Path) -> set[str]:
    return {line.strip() for line in filename.read_text(encoding="utf-8").splitlines() if line.strip()}


def parse_fasta_generator(filename: Path):
    current_seq_parts = []
    current_header = None
    with filename.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    yield current_header, "".join(current_seq_parts)
                current_header = line[1:]
                current_seq_parts = []
            else:
                current_seq_parts.append(line)
    if current_header is not None:
        yield current_header, "".join(current_seq_parts)


def main() -> None:
    fasta_input = INPUT_DIR / FASTA_INPUT_FILE
    names_file = INPUT_DIR / NAMES_FILE
    if not fasta_input.exists() or not names_file.exists():
        raise FileNotFoundError("输入 FASTA 或名称文件不存在。")
    names_to_extract = read_names_to_extract(names_file)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / FASTA_OUTPUT_FILE
    with output_file.open("w", encoding="utf-8") as out:
        for header, sequence in parse_fasta_generator(fasta_input):
            if header.split()[0] in names_to_extract:
                out.write(f">{header}\n{sequence}\n")
    print("按名称提取序列完成。")


if __name__ == "__main__":
    main()
