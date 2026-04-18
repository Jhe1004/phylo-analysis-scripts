#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
from pathlib import Path

from Bio import SeqIO


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/3_基因表格生成/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/3_基因表格生成/output"

GENBANK_SUBDIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/3_基因表格生成/input/gb_files"
# ==============================


OUTPUT_TSV_FILE = "gene_content_report.tsv"
OUTPUT_SUMMARY_FILE = "gene_summary_report.txt"


GENE_CLASSIFICATION = {
    "Self-replication": {
        "Ribosomal RNA genes": ["rrn23", "rrn16", "rrn5", "rrn4.5"],
        "Transfer RNA genes": [
            "trnA-UGC", "trnC-GCA", "trnD-GUC", "trnE-UUC", "trnF-GAA", "trnFm-CAU", "trnG-UCC",
            "trnG-GCC", "trnH-GUG", "trnI-CAU", "trnI-GAU", "trnK-UUU", "trnL-CAA", "trnL-UAA",
            "trnL-UAG", "trnM-CAU", "trnN-GUU", "trnP-UGG", "trnQ-UUG", "trnR-ACG", "trnR-UCU",
            "trnS-GCU", "trnS-GGA", "trnS-UGA", "trnT-GGU", "trnT-UGU", "trnV-GAC", "trnV-UAC",
            "trnW-CCA", "trnY-GUA",
        ],
        "Small subunit of ribosome": ["rps11", "rps12", "rps14", "rps15", "rps16", "rps18", "rps19", "rps2", "rps3", "rps4", "rps7", "rps8"],
        "Large subunit of ribosome": ["rpl14", "rpl16", "rpl2", "rpl20", "rpl22", "rpl23", "rpl32", "rpl33", "rpl36"],
        "DNA-dependent RNA polymerase": ["rpoA", "rpoB", "rpoC1", "rpoC2"],
        "Translational initiation factor": ["infA"],
    },
    "Genes for photosynthesis": {
        "Subunits of photosystem I": ["psaA", "psaB", "psaC", "psaI", "psaJ", "ycf4"],
        "Subunits of photosystem II": ["psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbT", "psbZ", "psbN"],
        "Subunits of cytochrome": ["petA", "petB", "petD", "petG", "petL", "petN"],
        "Subunits of ATP synthase": ["atpA", "atpB", "atpE", "atpF", "atpH", "atpI"],
        "Large subunit of Rubisco": ["rbcL"],
        "Subunits of NADH dehydrogenase": ["ndhA", "ndhB", "ndhC", "ndhD", "ndhE", "ndhF", "ndhG", "ndhH", "ndhI", "ndhJ", "ndhK"],
    },
    "Other genes": {
        "Maturase": ["matK"],
        "Envelope membrane protein": ["cemA"],
        "Subunit of acetyl-CoA": ["accD"],
        "C-type cytochrome synthesis gene": ["ccsA"],
        "Protease": ["clpP"],
        "Component of TIC complex": ["ycf1"],
        "Photosystem I assembly factor": ["ycf3"],
        "Essential unknown function": ["ycf2"],
    },
}


def main() -> None:
    input_dir = SCRIPT_DIR / INPUT_DIRECTORY / GENBANK_SUBDIRECTORY
    output_dir = SCRIPT_DIR / OUTPUT_DIRECTORY
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_classification_lower = {
        category: {group: [gene.lower() for gene in genes] for group, genes in groups.items()}
        for category, groups in GENE_CLASSIFICATION.items()
    }
    all_reference_genes = {
        gene
        for groups in gene_classification_lower.values()
        for genes in groups.values()
        for gene in genes
    }

    gb_files = [f for f in os.listdir(input_dir) if f.endswith((".gb", ".gbk", ".genbank"))]
    if not gb_files:
        raise FileNotFoundError(f"在 {input_dir} 中没有找到 GenBank 文件")

    presence_data = {}
    for filename in sorted(gb_files):
        filepath = input_dir / filename
        genes_in_file = set()
        for record in SeqIO.parse(str(filepath), "genbank"):
            for feature in record.features:
                if feature.type == "gene" and "gene" in feature.qualifiers:
                    genes_in_file.add(feature.qualifiers["gene"][0].strip().lower())
        presence_data[filename] = genes_in_file

    with open(output_dir / OUTPUT_TSV_FILE, "w", encoding="utf-8") as f_out:
        header = ["Category", "Gene group", "Gene name"] + sorted(presence_data.keys())
        f_out.write("\t".join(header) + "\n")
        for category, groups in GENE_CLASSIFICATION.items():
            for group_name, genes in groups.items():
                for gene in sorted(genes):
                    row_data = [category, group_name, gene]
                    gene_lower = gene.lower()
                    for filename in sorted(presence_data.keys()):
                        row_data.append("+" if gene_lower in presence_data[filename] else "-")
                    f_out.write("\t".join(row_data) + "\n")

    with open(output_dir / OUTPUT_SUMMARY_FILE, "w", encoding="utf-8") as f_summary:
        f_summary.write("Gene Content Analysis Summary\n")
        f_summary.write("=" * 30 + "\n\n")
        f_summary.write("--- Missing Genes (Present in Reference, Absent in Sample) ---\n")
        missing_summary = {
            gene: [Path(fn).stem for fn, found in presence_data.items() if gene not in found]
            for gene in sorted(all_reference_genes)
        }
        missing_summary = {gene: files for gene, files in missing_summary.items() if files}
        if not missing_summary:
            f_summary.write("No missing genes found across all samples based on the reference list.\n")
        else:
            for gene, files in missing_summary.items():
                if len(files) == len(gb_files):
                    f_summary.write(f"Gene: {gene:<10}\tMISSING in ALL {len(files)} samples.\n")
                else:
                    f_summary.write(f"Gene: {gene:<10}\tMissing in {len(files)} sample(s): {', '.join(files)}\n")
        f_summary.write("\n\n--- Extra Genes Found (Present in Sample, Absent in Reference) ---\n")
        all_found_genes = set.union(*presence_data.values())
        extra_genes = sorted(all_found_genes - all_reference_genes)
        if not extra_genes:
            f_summary.write("No extra genes were found that are not on the reference list.\n")
        else:
            for gene in extra_genes:
                found_in = [Path(fn).stem for fn, found in presence_data.items() if gene in found]
                f_summary.write(f"Gene: {gene:<15}\tFound in {len(found_in)} sample(s): {', '.join(found_in)}\n")

    print(f"完成。输出目录: {output_dir}")


if __name__ == "__main__":
    main()
