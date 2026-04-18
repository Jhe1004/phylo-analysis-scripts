#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import Counter
from pathlib import Path


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/1_统计各种数字/1_统计get_proteinortho结果中各个物种的基因数目/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/1_统计各种数字/1_统计get_proteinortho结果中各个物种的基因数目/output"
FILE_SUFFIX = "cds.fasta"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    fasta_files = sorted(INPUT_DIR.glob(f"*{FILE_SUFFIX}"))
    if not fasta_files:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到 *{FILE_SUFFIX} 文件。")
    total_species_counts = Counter()
    for filename in fasta_files:
        species_in_current_file = set()
        with filename.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    species_in_current_file.add(line.strip().split()[0][1:])
        total_species_counts.update(species_in_current_file)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    lines = ["species\tcount"]
    for species, count in total_species_counts.most_common():
        lines.append(f"{species}\t{count}")
    (OUTPUT_DIR / "species_counts.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("物种计数完成。")


if __name__ == "__main__":
    main()
