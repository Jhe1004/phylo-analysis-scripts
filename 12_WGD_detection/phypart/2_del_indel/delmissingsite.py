#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import os

from Bio import SeqIO


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/2_del_indel/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/2_del_indel/output"

MISSING_PROPORTION_THRESHOLD = 0.2
ALIGNMENT_CHUNK_SIZE = 50000
PROCESS_COUNT = max(1, min(8, os.cpu_count() or 1))
INPUT_EXTENSIONS = [".fasta", ".fas", ".fa"]
OUTPUT_EXTENSION = ".fas"
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_alignment(alignment_file: Path) -> list[tuple[str, str]]:
    records = []
    max_len = 0
    for record in SeqIO.parse(str(alignment_file), "fasta"):
        seq = str(record.seq)
        if not seq or all(char in "-?Nn" for char in seq):
            continue
        max_len = max(max_len, len(seq))
        records.append((record.id, seq))
    if not records:
        return []
    return [(seq_id, seq.ljust(max_len, "-")) for seq_id, seq in records]


def filter_alignment_columns(records: list[tuple[str, str]], threshold: float) -> list[tuple[str, str]]:
    if not records:
        return []
    seq_ids = [item[0] for item in records]
    sequences = [item[1] for item in records]
    seq_count = len(sequences)
    kept_columns = []
    for column in zip(*sequences):
        gap_count = sum(1 for char in column if char not in {"A", "a", "T", "t", "C", "c", "G", "g"})
        if gap_count / seq_count < threshold:
            kept_columns.append(column)
    if not kept_columns:
        return []
    filtered_sequences = ["".join(chars) for chars in zip(*kept_columns)]
    return [
        (seq_id, filtered_seq)
        for seq_id, filtered_seq in zip(seq_ids, filtered_sequences)
        if filtered_seq and not all(char in "-?" for char in filtered_seq)
    ]


def split_records(records: list[tuple[str, str]], chunk_size: int) -> list[list[tuple[str, str]]]:
    if not records:
        return []
    length = len(records[0][1])
    if length <= chunk_size:
        return [records]
    chunks = []
    for left in range(0, length, chunk_size):
        right = min(length, left + chunk_size)
        chunks.append([(seq_id, seq[left:right]) for seq_id, seq in records])
    return chunks


def merge_chunks(chunks: list[list[tuple[str, str]]]) -> list[tuple[str, str]]:
    merged: dict[str, list[str]] = {}
    for chunk in chunks:
        for seq_id, seq in chunk:
            merged.setdefault(seq_id, []).append(seq)
    return [(seq_id, "".join(parts)) for seq_id, parts in merged.items()]


def write_alignment(records: list[tuple[str, str]], output_file: Path) -> None:
    with open(output_file, "w", encoding="utf-8") as handle:
        for seq_id, seq in records:
            handle.write(f">{seq_id}\n{seq}\n")


def process_alignment(task: tuple[Path, Path, float, int]) -> tuple[str, int, int]:
    input_file, output_file, threshold, chunk_size = task
    records = read_alignment(input_file)
    if not records:
        return input_file.name, 0, 0
    chunks = split_records(records, chunk_size)
    filtered_chunks = [filter_alignment_columns(chunk, threshold) for chunk in chunks]
    filtered_chunks = [chunk for chunk in filtered_chunks if chunk]
    if not filtered_chunks:
        return input_file.name, len(records[0][1]), 0
    merged = merge_chunks(filtered_chunks)
    write_alignment(merged, output_file)
    return input_file.name, len(records[0][1]), len(merged[0][1]) if merged else 0


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    ensure_directory(output_dir)

    alignment_files = []
    for extension in INPUT_EXTENSIONS:
        alignment_files.extend(sorted(input_dir.glob(f"*{extension}")))
    if not alignment_files:
        raise FileNotFoundError(f"未在 {input_dir} 中找到待处理比对文件")

    tasks = [
        (
            alignment_file,
            output_dir / f"{alignment_file.stem}{OUTPUT_EXTENSION}",
            MISSING_PROPORTION_THRESHOLD,
            ALIGNMENT_CHUNK_SIZE,
        )
        for alignment_file in alignment_files
    ]

    if PROCESS_COUNT > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=PROCESS_COUNT) as executor:
            results = list(executor.map(process_alignment, tasks))
    else:
        results = [process_alignment(task) for task in tasks]

    summary_file = output_dir / "delmissing_summary.tsv"
    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write("file\toriginal_length\tfiltered_length\n")
        for file_name, original_length, filtered_length in results:
            handle.write(f"{file_name}\t{original_length}\t{filtered_length}\n")

    print(f"完成矩阵数: {len(results)}")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
