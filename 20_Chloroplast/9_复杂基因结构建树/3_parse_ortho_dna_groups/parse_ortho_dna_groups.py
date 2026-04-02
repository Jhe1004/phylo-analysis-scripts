#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
TSV_GLOB_PATTERN = "proteinortho_*.tsv"
FASTA_SUBDIRECTORY = "source_fastas"
SINGLE_COPY_ONLY = True
MIN_SPECIES_PER_OG = 30


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY
SOURCE_FASTA_DIR = INPUT_DIR / FASTA_SUBDIRECTORY


def load_all_sequences(species_files: list[str]) -> dict[str, SeqRecord]:
    all_seq_dict: dict[str, SeqRecord] = {}
    required_files = set()
    for file_name in species_files:
        if file_name.endswith(".fasta"):
            required_files.add(file_name)
        elif file_name.endswith(".gb") or file_name.endswith(".gbk"):
            required_files.add(f"{Path(file_name).stem}.fasta")
    for fasta_name in required_files:
        fasta_path = SOURCE_FASTA_DIR / fasta_name
        if not fasta_path.exists():
            print(f"警告: 未找到源 FASTA 文件 {fasta_path}")
            continue
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            all_seq_dict[record.id] = record
    return all_seq_dict


def process_one_tsv(tsv_file: Path, output_dir: Path, all_sequences_db: dict[str, SeqRecord]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    species_files_in_header: list[str] = []
    with tsv_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("# Species"):
                species_files_in_header = line.strip().split("\t")[3:]
                break
    if not species_files_in_header:
        raise RuntimeError(f"无法从 {tsv_file} 表头解析物种列表。")

    og_counter = 0
    og_written = 0
    with tsv_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            try:
                num_species = int(cols[0])
                num_genes = int(cols[1])
            except (ValueError, IndexError):
                continue
            if num_species < MIN_SPECIES_PER_OG:
                continue
            if SINGLE_COPY_ONLY and num_species != num_genes:
                continue

            og_counter += 1
            og_name = f"OG_{og_counter:05d}"
            og_records: list[SeqRecord] = []
            for index, species_file in enumerate(species_files_in_header):
                species_id = Path(species_file).stem
                cell_data = cols[index + 3]
                if cell_data == "*":
                    continue
                seq_ids = cell_data.split(",")
                if SINGLE_COPY_ONLY and len(seq_ids) > 1:
                    continue
                for part_index, seq_id in enumerate(seq_ids, start=1):
                    if seq_id not in all_sequences_db:
                        continue
                    original = all_sequences_db[seq_id]
                    record_id = species_id if len(seq_ids) == 1 else f"{species_id}_p{part_index}"
                    og_records.append(SeqRecord(original.seq, id=record_id, description=f"original_id={seq_id}"))

            if og_records:
                SeqIO.write(og_records, str(output_dir / f"{og_name}.fasta"), "fasta")
                og_written += 1
    print(f"{tsv_file.name}: 共写出 {og_written} 个同源组。")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tsv_files = sorted(INPUT_DIR.glob(TSV_GLOB_PATTERN))
    if not tsv_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到匹配 {TSV_GLOB_PATTERN} 的文件。")

    species_files_in_header: list[str] = []
    with tsv_files[0].open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("# Species"):
                species_files_in_header = line.strip().split("\t")[3:]
                break
    if not species_files_in_header:
        raise RuntimeError(f"无法从 {tsv_files[0]} 解析物种列表。")

    all_sequences_db = load_all_sequences(species_files_in_header)
    if not all_sequences_db:
        raise RuntimeError("未能加载任何源序列。")

    for tsv_file in tsv_files:
        suffix = tsv_file.stem.replace("proteinortho_", "")
        process_one_tsv(tsv_file, OUTPUT_DIR / f"{suffix}_groups", all_sequences_db)
    print("Proteinortho 同源组解析完成。")


if __name__ == "__main__":
    main()
