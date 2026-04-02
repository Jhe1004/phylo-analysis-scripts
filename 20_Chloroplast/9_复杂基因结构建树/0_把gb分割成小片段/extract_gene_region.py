#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
GENBANK_SUBDIRECTORY = "gb_files"
FEATURES_TO_EXTRACT = ["CDS", "rRNA", "tRNA"]


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY / GENBANK_SUBDIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def split_genome_to_granular_segments() -> None:
    gb_files = sorted(INPUT_DIR.glob("*.gb")) + sorted(INPUT_DIR.glob("*.gbk"))
    if not gb_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 .gb 或 .gbk 文件。")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for gb_file in gb_files:
        base_name = gb_file.stem
        output_filename = OUTPUT_DIR / f"{base_name}.fasta"
        record = SeqIO.read(str(gb_file), "genbank")
        fasta_segments = []
        current_position = 0
        total_length = len(record.seq)
        all_parts = []

        for feature in record.features:
            if feature.type not in FEATURES_TO_EXTRACT:
                continue
            gene_name = feature.qualifiers.get("gene", [feature.type])[0]
            if len(feature.location.parts) > 1:
                for index, part in enumerate(feature.location.parts, start=1):
                    part_name = f"{gene_name}_exon{index}"
                    all_parts.append((int(part.start), int(part.end), part_name, feature.type, part))
            else:
                part = feature.location
                all_parts.append((int(part.start), int(part.end), gene_name, feature.type, part))

        sorted_parts = sorted(all_parts, key=lambda item: item[0])
        for part_start, part_end, part_name, part_type, part_location in sorted_parts:
            if part_start > current_position:
                gap_seq = record.seq[current_position:part_start]
                gap_id = f"{base_name}_IGS_{current_position + 1}-{part_start}"
                fasta_segments.append(SeqRecord(gap_seq, id=gap_id, description=""))
            part_seq = part_location.extract(record.seq)
            feature_id = f"{base_name}_{part_name}_{part_start + 1}-{part_end}"
            fasta_segments.append(SeqRecord(part_seq, id=feature_id, description=""))
            current_position = max(current_position, part_end)

        if current_position < total_length:
            last_gap_seq = record.seq[current_position:total_length]
            last_gap_id = f"{base_name}_IGS_{current_position + 1}-{total_length}"
            fasta_segments.append(SeqRecord(last_gap_seq, id=last_gap_id, description=""))

        SeqIO.write(fasta_segments, str(output_filename), "fasta")
        print(f"已生成: {output_filename}")


if __name__ == "__main__":
    split_genome_to_granular_segments()
