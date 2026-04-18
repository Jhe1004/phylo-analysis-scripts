#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord


# =========================
# 用户配置区
# =========================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/8_叶绿体用cds建树/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/8_叶绿体用cds建树/output"
CONDA_ENV_NAME = "trinity_env"

GENBANK_SUBDIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/8_叶绿体用cds建树/input/gb_files"
REFERENCE_FILENAME = ""
E_VALUE_THRESHOLD = 1e-10
MIN_OUTPUT_SEQUENCES = 2
KEEP_TEMPORARY_FILES = False


SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY
DEPENDENCIES_DIR = SCRIPT_DIR / "dependencies"


def resolve_user_path(path_text: str) -> Path:
    return SCRIPT_DIR / path_text


def resolve_conda_executable(env_name: str, executable_name: str) -> Path | None:
    candidate_dirs: list[Path] = []
    conda_prefix = Path.home() / "miniconda3" / "envs" / env_name / "bin"
    mamba_prefix = Path.home() / "mambaforge" / "envs" / env_name / "bin"
    anaconda_prefix = Path.home() / "anaconda3" / "envs" / env_name / "bin"
    candidate_dirs.extend([conda_prefix, mamba_prefix, anaconda_prefix])
    for directory in candidate_dirs:
        candidate = directory / executable_name
        if candidate.exists():
            return candidate
    return None


def find_executable(executable_name: str) -> Path:
    local_candidate = DEPENDENCIES_DIR / "bin" / executable_name
    if local_candidate.exists():
        return local_candidate
    conda_candidate = resolve_conda_executable(CONDA_ENV_NAME, executable_name)
    if conda_candidate is not None:
        return conda_candidate
    system_candidate = shutil.which(executable_name)
    if system_candidate:
        return Path(system_candidate)
    raise FileNotFoundError(
        f"未找到可执行文件 '{executable_name}'。查找顺序: dependencies/bin -> conda 环境 '{CONDA_ENV_NAME}' -> PATH"
    )


def sanitize_filename(name: str) -> str:
    return "".join(c for c in name if c.isalnum() or c in ("_", "-", ".")).rstrip()


def sanitize_sample_name(name: str) -> str:
    safe = "".join(c if c.isalnum() or c in ("_", "-") else "_" for c in name)
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe.strip("_")


def load_genbank_files(input_dir: Path) -> list[Path]:
    gb_dir = input_dir / GENBANK_SUBDIRECTORY
    gb_files = sorted(gb_dir.glob("*.gb")) + sorted(gb_dir.glob("*.gbk"))
    if not gb_files:
        raise FileNotFoundError(
            f"未在 {gb_dir} 中找到 .gb 或 .gbk 文件。请把 GenBank 文件放到 {INPUT_DIRECTORY}/{GENBANK_SUBDIRECTORY}/"
        )
    return gb_files


def choose_reference(gb_files: list[Path]) -> Path:
    if REFERENCE_FILENAME:
        reference_path = INPUT_DIR / GENBANK_SUBDIRECTORY / REFERENCE_FILENAME
        if not reference_path.exists():
            raise FileNotFoundError(f"指定的参考文件不存在: {reference_path}")
        return reference_path
    return gb_files[0]


def create_blast_database(genbank_file: Path, fasta_path: Path, db_prefix: Path, makeblastdb: Path) -> None:
    SeqIO.convert(str(genbank_file), "genbank", str(fasta_path), "fasta")
    command = [
        str(makeblastdb),
        "-in",
        str(fasta_path),
        "-dbtype",
        "nucl",
        "-out",
        str(db_prefix),
    ]
    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def get_best_feature(record, gene_name: str, candidates: Iterable | None = None):
    if candidates is None:
        candidates = [
            feature
            for feature in record.features
            if "gene" in feature.qualifiers and feature.qualifiers["gene"][0] == gene_name
        ]
    valid = [
        feature
        for feature in candidates
        if "pseudo" not in feature.qualifiers and "pseudogene" not in feature.qualifiers
    ]
    if not valid:
        return None, None
    priority_order = ["CDS", "rRNA", "tRNA", "gene"]

    def sort_key(feature):
        try:
            priority = priority_order.index(feature.type)
        except ValueError:
            priority = len(priority_order)
        return (priority, -len(feature))

    valid = sorted(valid, key=sort_key)
    best = valid[0]
    return best, best.type


def find_ortholog_feature_in_hit(target_record, gene_name: str, hit_start: int, hit_end: int):
    hit_location = set(range(hit_start, hit_end + 1))
    candidates = []
    for feature in target_record.features:
        if "gene" not in feature.qualifiers:
            continue
        if feature.qualifiers["gene"][0] != gene_name:
            continue
        feature_location = set(range(int(feature.location.start), int(feature.location.end)))
        if hit_location.intersection(feature_location):
            candidates.append(feature)
    if not candidates:
        return None
    best_feature, _ = get_best_feature(target_record, gene_name, candidates=candidates)
    if best_feature is None:
        return None
    return best_feature.extract(target_record.seq)


def main() -> None:
    makeblastdb = find_executable("makeblastdb")
    blastn = find_executable("blastn")
    print(f"使用 makeblastdb: {makeblastdb}")
    print(f"使用 blastn: {blastn}")

    gb_files = load_genbank_files(INPUT_DIR)
    reference_file = choose_reference(gb_files)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ortholog_dir = OUTPUT_DIR / "ortholog_fastas"
    ortholog_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=str(OUTPUT_DIR), prefix="tmp_blast_") as temp_name:
        temp_dir = Path(temp_name)
        db_paths: dict[str, Path] = {}
        all_records = {}

        print("开始创建 BLAST 数据库并读取 GenBank 记录...")
        for gb_file in gb_files:
            sample_name = sanitize_sample_name(gb_file.stem)
            fasta_path = temp_dir / f"{sample_name}.fasta"
            db_prefix = temp_dir / sample_name
            create_blast_database(gb_file, fasta_path, db_prefix, makeblastdb)
            db_paths[sample_name] = db_prefix
            all_records[sample_name] = SeqIO.read(str(gb_file), "genbank")

        reference_name = sanitize_sample_name(reference_file.stem)
        reference_record = all_records[reference_name]
        all_gene_names = sorted(
            {
                feature.qualifiers["gene"][0]
                for feature in reference_record.features
                if "gene" in feature.qualifiers
            }
        )

        summary_lines = [
            f"输入目录: {INPUT_DIR / GENBANK_SUBDIRECTORY}",
            f"输出目录: {ortholog_dir}",
            f"参考文件: {reference_file.name}",
            f"E-value 阈值: {E_VALUE_THRESHOLD}",
            "",
        ]

        for gene_name_raw in all_gene_names:
            best_feature, extraction_type = get_best_feature(reference_record, gene_name_raw)
            if best_feature is None:
                continue

            gene_name = sanitize_filename(gene_name_raw)
            if not gene_name:
                continue

            ref_gene_seq = best_feature.extract(reference_record.seq)
            query_fasta = temp_dir / f"query_{gene_name}.fasta"
            with query_fasta.open("w", encoding="utf-8") as handle:
                handle.write(f">{gene_name}_ref\n{ref_gene_seq}\n")

            orthologs = [
                SeqRecord(
                    ref_gene_seq,
                    id=reference_name,
                    description=f"gene={gene_name_raw} type={extraction_type}",
                )
            ]

            for sample_name, db_prefix in db_paths.items():
                if sample_name == reference_name:
                    continue
                target_record = all_records[sample_name]
                blast_xml = temp_dir / f"{gene_name}_{sample_name}.xml"
                blastn_cline = NcbiblastnCommandline(
                    cmd=str(blastn),
                    query=str(query_fasta),
                    db=str(db_prefix),
                    evalue=E_VALUE_THRESHOLD,
                    outfmt=5,
                    out=str(blast_xml),
                )
                blastn_cline()

                try:
                    with blast_xml.open() as result_handle:
                        blast_record = next(NCBIXML.parse(result_handle), None)
                    if not blast_record or not blast_record.alignments:
                        continue
                    hsp = blast_record.alignments[0].hsps[0]
                    if hsp.expect > E_VALUE_THRESHOLD:
                        continue
                    start = min(hsp.sbjct_start, hsp.sbjct_end)
                    end = max(hsp.sbjct_start, hsp.sbjct_end)
                    extracted = find_ortholog_feature_in_hit(target_record, gene_name_raw, start, end)
                    if extracted is None:
                        continue
                    orthologs.append(
                        SeqRecord(
                            extracted,
                            id=sample_name,
                            description=f"ortholog_of={gene_name_raw}",
                        )
                    )
                except Exception:
                    continue

            if len(orthologs) >= MIN_OUTPUT_SEQUENCES:
                output_fasta = ortholog_dir / f"{gene_name}.fasta"
                SeqIO.write(orthologs, str(output_fasta), "fasta")
                summary_lines.append(f"{gene_name}\t{len(orthologs)}")

        (OUTPUT_DIR / "ortholog_summary.tsv").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

        if KEEP_TEMPORARY_FILES:
            kept_dir = OUTPUT_DIR / "kept_temp_blast"
            if kept_dir.exists():
                shutil.rmtree(kept_dir)
            shutil.copytree(temp_dir, kept_dir)

    print("直系同源基因提取完成。")


if __name__ == "__main__":
    main()
