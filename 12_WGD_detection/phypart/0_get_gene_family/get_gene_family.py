#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import csv

from Bio import SeqIO


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

PROTEINORTHO_FILE = "myproject.proteinortho.tsv"
SEQUENCE_SUBDIRECTORY = "sequences"
MIN_SPECIES_COUNT = 4
SKIP_SINGLE_COPY_FAMILIES = False
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_sequence_dicts(sequence_dir: Path) -> tuple[dict[str, dict], dict[str, dict]]:
    pep_dicts: dict[str, dict] = {}
    cds_dicts: dict[str, dict] = {}
    for sequence_file in sorted(sequence_dir.iterdir()):
        if not sequence_file.is_file():
            continue
        suffix = sequence_file.name
        if suffix.endswith(".pep"):
            pep_dicts[sequence_file.name] = SeqIO.to_dict(SeqIO.parse(str(sequence_file), "fasta"))
        elif suffix.endswith(".cds"):
            cds_dicts[sequence_file.name] = SeqIO.to_dict(SeqIO.parse(str(sequence_file), "fasta"))
    return pep_dicts, cds_dicts


def get_partner_file_name(species_file_name: str, target_suffix: str) -> str:
    if species_file_name.endswith(".pep"):
        return species_file_name[:-4] + target_suffix
    if species_file_name.endswith(".cds"):
        return species_file_name[:-4] + target_suffix
    raise ValueError(f"无法识别物种序列文件后缀: {species_file_name}")


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    sequence_dir = input_dir / SEQUENCE_SUBDIRECTORY
    family_dir = output_dir / "gene_families"
    ensure_directory(family_dir)

    proteinortho_file = input_dir / PROTEINORTHO_FILE
    if not proteinortho_file.exists():
        raise FileNotFoundError(f"未找到 proteinortho 结果文件: {proteinortho_file}")

    pep_dicts, cds_dicts = load_sequence_dicts(sequence_dir)
    if not pep_dicts or not cds_dicts:
        raise FileNotFoundError(f"未在 {sequence_dir} 中找到 .pep/.cds 文件")

    family_count = 0
    exported_count = 0
    summary_rows: list[list[str]] = []

    with open(proteinortho_file, "r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        species_columns = header[3:]
        for row_index, row in enumerate(reader, start=1):
            family_count += 1
            reference_count = int(float(row[0]))
            species_count = int(float(row[1]))
            if SKIP_SINGLE_COPY_FAMILIES and reference_count == species_count:
                continue
            if species_count < MIN_SPECIES_COUNT:
                continue

            family_name = f"ortho{row_index}"
            pep_output = family_dir / f"{family_name}_pep.fasta"
            cds_output = family_dir / f"{family_name}_cds.fasta"
            sequence_written = 0

            with open(pep_output, "w", encoding="utf-8") as pep_handle, open(
                cds_output, "w", encoding="utf-8"
            ) as cds_handle:
                for column_index, species_file_name in enumerate(species_columns, start=3):
                    cell = row[column_index]
                    if cell == "*":
                        continue
                    pep_file_name = get_partner_file_name(species_file_name, ".pep")
                    cds_file_name = get_partner_file_name(species_file_name, ".cds")
                    pep_records = pep_dicts.get(pep_file_name)
                    cds_records = cds_dicts.get(cds_file_name)
                    if pep_records is None or cds_records is None:
                        raise FileNotFoundError(
                            f"缺少配套序列文件: {species_file_name} -> {pep_file_name} / {cds_file_name}"
                        )
                    for sequence_id in cell.split(","):
                        pep_record = pep_records.get(sequence_id)
                        cds_record = cds_records.get(sequence_id)
                        if pep_record is None or cds_record is None:
                            raise KeyError(f"序列 ID 不存在: {species_file_name} ++ {sequence_id}")
                        header_id = f"{species_file_name}++{sequence_id}"
                        pep_handle.write(f">{header_id}\n{str(pep_record.seq)}\n")
                        cds_handle.write(f">{header_id}\n{str(cds_record.seq)}\n")
                        sequence_written += 1

            if sequence_written == 0:
                pep_output.unlink(missing_ok=True)
                cds_output.unlink(missing_ok=True)
                continue

            exported_count += 1
            summary_rows.append([family_name, str(species_count), str(sequence_written)])

    summary_file = output_dir / "gene_family_summary.tsv"
    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write("family\tspecies_count\tsequence_count\n")
        for row in summary_rows:
            handle.write("\t".join(row) + "\n")

    print(f"读取基因家族数: {family_count}")
    print(f"导出基因家族数: {exported_count}")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
