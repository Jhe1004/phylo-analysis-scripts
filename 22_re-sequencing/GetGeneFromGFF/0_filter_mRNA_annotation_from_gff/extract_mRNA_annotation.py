#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/22_re-sequencing/GetGeneFromGFF/0_filter_mRNA_annotation_from_gff/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/22_re-sequencing/GetGeneFromGFF/0_filter_mRNA_annotation_from_gff/output"
INPUT_GFF_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/22_re-sequencing/GetGeneFromGFF/0_filter_mRNA_annotation_from_gff/input/ref.gff"
OUTPUT_GFF_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/22_re-sequencing/GetGeneFromGFF/0_filter_mRNA_annotation_from_gff/output/filtered_mrna.gff"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def read_gff(filename: Path) -> list[list[str]]:
    mrna_annotations = []
    with filename.open("r", encoding="utf-8") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 9 and fields[2] == "mRNA":
                mrna_annotations.append(fields)
    return mrna_annotations


def remove_overlapping_annotations(annotations: list[list[str]]) -> list[list[str]]:
    non_overlapping = []
    sorted_annotations = sorted(annotations, key=lambda item: (item[0], int(item[3]), int(item[4])))
    current_chr = ""
    current_end = 0
    for annotation in sorted_annotations:
        chromosome, start, end = annotation[0], int(annotation[3]), int(annotation[4])
        if chromosome != current_chr or start > current_end:
            non_overlapping.append(annotation)
            current_chr = chromosome
            current_end = end
        else:
            current_end = max(current_end, end)
    return non_overlapping


def main() -> None:
    gff_file = INPUT_DIR / INPUT_GFF_FILE
    if not gff_file.exists():
        raise FileNotFoundError(f"未找到输入 GFF 文件: {gff_file}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / OUTPUT_GFF_FILE

    mrna_annotations = read_gff(gff_file)
    filtered_annotations = remove_overlapping_annotations(mrna_annotations)

    with output_file.open("w", encoding="utf-8") as handle:
        for annotation in filtered_annotations:
            handle.write("\t".join(annotation) + "\n")

    print(f"已输出过滤后的 mRNA 注释: {output_file}")


if __name__ == "__main__":
    main()
