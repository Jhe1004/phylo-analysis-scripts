#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import functools
import multiprocessing
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/21_annotation/annotate_fasta_genes/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/21_annotation/annotate_fasta_genes/output"
CONDA_ENV_NAME = "trinity_env"

PEP_SUBDIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/21_annotation/annotate_fasta_genes/input/raw"
REFERENCE_FASTA_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/21_annotation/annotate_fasta_genes/input/reference/tair_filter.fasta.transdecoder.pep"

PEP_FILENAME_PATTERN = "*_pep_maffted.fas"
PEP_SUFFIX = "_pep_maffted.fas"
CDS_FILENAME_TEMPLATE = "{}_cds_maffted.fas"

DB_NAME = "arabidopsis_db"
EVALUE_THRESHOLD = 1e-5
MAX_TARGET_SEQS = 1
PROCESS_COUNT = 0
THREADS_PER_BLAST = 1


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY
DEPENDENCIES_DIR = SCRIPT_DIR / "dependencies"


def resolve_conda_executable(env_name: str, executable_name: str) -> Path | None:
    for prefix in (
        Path.home() / "miniconda3" / "envs" / env_name / "bin",
        Path.home() / "mambaforge" / "envs" / env_name / "bin",
        Path.home() / "anaconda3" / "envs" / env_name / "bin",
    ):
        candidate = prefix / executable_name
        if candidate.exists():
            return candidate
    return None


def find_executable(name: str) -> Path:
    local_candidate = DEPENDENCIES_DIR / "bin" / name
    if local_candidate.exists():
        return local_candidate
    conda_candidate = resolve_conda_executable(CONDA_ENV_NAME, name)
    if conda_candidate is not None:
        return conda_candidate
    system_candidate = shutil.which(name)
    if system_candidate:
        return Path(system_candidate)
    raise FileNotFoundError(
        f"未找到可执行文件 '{name}'。查找顺序: dependencies/bin -> conda 环境 '{CONDA_ENV_NAME}' -> PATH"
    )


def make_blast_db(arabidopsis_fasta: Path, db_name: Path, makeblastdb: Path) -> bool:
    db_files = [db_name.with_suffix(f".{ext}") for ext in ["phr", "pin", "psq"]]
    if all(path.exists() for path in db_files):
        print(f"BLAST 数据库 '{db_name}' 已存在，跳过创建。")
        return True
    cmd = [str(makeblastdb), "-in", str(arabidopsis_fasta), "-dbtype", "prot", "-out", str(db_name)]
    subprocess.run(cmd, check=True, capture_output=True, text=True, encoding="utf-8")
    return True


def run_blastp_in_tempdir(
    cleaned_query_fasta: Path,
    db_name: Path,
    temp_blast_output: Path,
    blastp: Path,
    evalue: float,
    max_target_seqs: int,
    num_threads: int,
) -> bool:
    cmd = [
        str(blastp),
        "-query",
        str(cleaned_query_fasta),
        "-db",
        str(db_name),
        "-out",
        str(temp_blast_output),
        "-evalue",
        str(evalue),
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(num_threads),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True, encoding="utf-8")
    return True


def parse_blast_output(blast_output: Path):
    columns = [
        "query_id",
        "subject_id",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
    if not blast_output.exists() or blast_output.stat().st_size == 0:
        return None
    try:
        blast_df = pd.read_csv(blast_output, sep="\t", names=columns)
    except pd.errors.EmptyDataError:
        return None
    return None if blast_df.empty else blast_df


def determine_best_gene(blast_df) -> str | None:
    if blast_df is None or blast_df.empty:
        return None
    gene_scores = blast_df.groupby("subject_id")["bitscore"].sum().reset_index()
    if gene_scores.empty:
        return None
    best_gene = gene_scores.sort_values(by="bitscore", ascending=False).iloc[0]["subject_id"]
    return best_gene.split(".")[0] if "." in best_gene else best_gene


def find_corresponding_cds_file(pep_fasta_path: Path, input_dir: Path) -> Path | None:
    base_name = pep_fasta_path.name
    if not base_name.endswith(PEP_SUFFIX):
        return None
    stem = base_name[: -len(PEP_SUFFIX)]
    if not stem:
        return None
    expected_cds_filename = CDS_FILENAME_TEMPLATE.format(stem)
    cds_search_path = input_dir / expected_cds_filename
    return cds_search_path if cds_search_path.is_file() else None


def clean_fasta_for_blast(original_fasta_path: Path, cleaned_fasta_path: Path) -> bool:
    cleaned_records = []
    for record in SeqIO.parse(str(original_fasta_path), "fasta"):
        cleaned_seq_string = str(record.seq).replace("-", "").replace("*", "")
        record.seq = type(record.seq)(cleaned_seq_string)
        cleaned_records.append(record)
    with cleaned_fasta_path.open("w", encoding="utf-8") as outfile:
        SeqIO.write(cleaned_records, outfile, "fasta")
    return True


def process_single_fasta(
    original_pep_fasta_path: Path,
    input_dir: Path,
    output_dir: Path,
    db_name: Path,
    blastp: Path,
    evalue: float,
    max_target_seqs: int,
    num_threads_per_blast: int,
):
    process_id = multiprocessing.current_process().pid
    base_name = original_pep_fasta_path.name
    temp_dir = None
    best_gene_found = None
    status = "失败"
    error_message = ""

    try:
        temp_dir = Path(tempfile.mkdtemp(prefix=f"blast_{process_id}_{base_name}_", dir=str(output_dir.parent)))
        cleaned_fasta_path = temp_dir / f"cleaned_{base_name}"
        clean_fasta_for_blast(original_pep_fasta_path, cleaned_fasta_path)

        temp_blast_output = temp_dir / f"{base_name}.blastp.out"
        run_blastp_in_tempdir(
            cleaned_fasta_path,
            db_name,
            temp_blast_output,
            blastp,
            evalue,
            max_target_seqs,
            num_threads_per_blast,
        )

        blast_df = parse_blast_output(temp_blast_output)
        if blast_df is None:
            return base_name, "无命中", None, "无 BLAST 命中或解析错误"

        best_gene = determine_best_gene(blast_df)
        if not best_gene:
            return base_name, "无最佳基因", None, "无法确定最佳基因"

        best_gene_found = best_gene
        pep_output_path = output_dir / f"{best_gene}_pep{original_pep_fasta_path.suffix}"
        shutil.copy2(original_pep_fasta_path, pep_output_path)

        original_cds_path = find_corresponding_cds_file(original_pep_fasta_path, input_dir)
        if original_cds_path:
            cds_output_path = output_dir / f"{best_gene}_cds{original_cds_path.suffix}"
            shutil.copy2(original_cds_path, cds_output_path)
        else:
            error_message = "未找到对应 CDS 文件"

        status = "成功"
    except Exception as exc:
        status = "失败"
        error_message = str(exc)
    finally:
        if temp_dir and temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)
    return base_name, status, best_gene_found, error_message


def run_pipeline() -> None:
    input_dir = INPUT_DIR / PEP_SUBDIRECTORY
    reference_fasta = INPUT_DIR / REFERENCE_FASTA_FILE
    output_dir = OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.is_dir():
        raise FileNotFoundError(f"输入目录无效: {input_dir}")
    if not reference_fasta.is_file():
        raise FileNotFoundError(f"参考蛋白 FASTA 不存在: {reference_fasta}")

    makeblastdb = find_executable("makeblastdb")
    blastp = find_executable("blastp")
    print(f"使用 makeblastdb: {makeblastdb}")
    print(f"使用 blastp: {blastp}")

    db_name = output_dir / DB_NAME
    make_blast_db(reference_fasta, db_name, makeblastdb)

    pep_fasta_files = sorted(input_dir.glob(PEP_FILENAME_PATTERN))
    if not pep_fasta_files:
        raise FileNotFoundError(f"未找到匹配 {PEP_FILENAME_PATTERN} 的蛋白 FASTA 文件。")

    if PROCESS_COUNT == 0:
        cpu_count = multiprocessing.cpu_count()
        num_workers = cpu_count // 2 if cpu_count > 1 else 1
    else:
        num_workers = PROCESS_COUNT
    num_workers = max(1, min(num_workers, len(pep_fasta_files)))

    start_time = time.time()
    worker_func = functools.partial(
        process_single_fasta,
        input_dir=input_dir,
        output_dir=output_dir,
        db_name=db_name,
        blastp=blastp,
        evalue=EVALUE_THRESHOLD,
        max_target_seqs=MAX_TARGET_SEQS,
        num_threads_per_blast=THREADS_PER_BLAST,
    )

    if num_workers > 1:
        with multiprocessing.Pool(processes=num_workers) as pool:
            results = pool.map(worker_func, pep_fasta_files)
    else:
        results = [worker_func(path) for path in pep_fasta_files]

    end_time = time.time()
    summary_lines = ["file\tstatus\tbest_gene\terror_message"]
    for base_name, status, best_gene, error_message in results:
        summary_lines.append(f"{base_name}\t{status}\t{best_gene or ''}\t{error_message}")
    summary_lines.append(f"\n总耗时\t{end_time - start_time:.2f} 秒")
    (output_dir / "annotation_summary.tsv").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("FASTA 注释完成。")


if __name__ == "__main__":
    run_pipeline()
