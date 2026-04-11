#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import glob
import logging
import os
import sys
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def ensure_exists(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"未找到{description}: {path}")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


# =============================
# 参数配置区
# =============================
CONDA_ENV_NAME = "reseq"
REFERENCE_GENOME = "input/reference/ref.fasta"
VCF_DIRECTORY = "output/vcf_output"
OUTPUT_DIRECTORY = "output/consensus_fasta_output"

PARALLEL_SAMPLES = 10
MIN_GQ = 20
MIN_DP = 3
USE_IUPAC = True


def setup_logger(script_name: str) -> logging.Logger:
    log_file = SCRIPT_DIR / f"{script_name}.log"
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    logger.info(f"日志文件: {log_file}")
    return logger


LOGGER = setup_logger("04_vcf_to_fasta")

try:
    import pysam
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError as exc:
    raise ImportError("需要 pysam 和 biopython") from exc

IUPAC_CODES = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A", "C", "G", "T"]): "N",
}


def get_iupac_code(alleles: list[str]) -> str:
    alleles_set = frozenset([a.upper() for a in alleles if a.upper() in {"A", "C", "G", "T"}])
    if len(alleles_set) == 0:
        return "N"
    if len(alleles_set) == 1:
        return list(alleles_set)[0]
    return IUPAC_CODES.get(alleles_set, "N")


def parse_reference_fasta(ref_path: str) -> OrderedDict[str, list[str]]:
    sequences: OrderedDict[str, list[str]] = OrderedDict()
    for record in SeqIO.parse(ref_path, "fasta"):
        sequences[record.id] = list(str(record.seq).upper())
    return sequences


def vcf_to_consensus(
    vcf_path: str,
    reference_seqs: OrderedDict[str, list[str]],
    output_folder: str,
    use_iupac: bool,
    min_gq: int,
    min_dp: int,
) -> list[str]:
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    sample_seqs = {
        sample: {contig: ["N"] * len(seq) for contig, seq in reference_seqs.items()}
        for sample in samples
    }

    for record in vcf:
        contig = record.chrom
        start = record.start
        stop = record.stop
        if contig not in reference_seqs:
            continue
        ref_allele = record.ref
        alt_alleles = record.alts if record.alts else []

        for sample in samples:
            gt = record.samples[sample].get("GT")
            if gt is None or None in gt:
                continue
            gq = record.samples[sample].get("GQ", 99)
            dp = record.samples[sample].get("DP", 100)
            if (gq is not None and gq < min_gq) or (dp is not None and dp < min_dp):
                continue

            alleles_list = [ref_allele] + list(alt_alleles)
            try:
                called_alleles = [alleles_list[i] for i in gt]
            except (IndexError, TypeError):
                continue

            called_alleles_clean = [a for a in called_alleles if all(c in "ACGTacgt" for c in a)]
            if not called_alleles_clean:
                if any(a == "<NON_REF>" for a in called_alleles) and len(set(gt)) == 1 and gt[0] == 0:
                    called_alleles_clean = [ref_allele]
                else:
                    continue

            unique_alleles = set(called_alleles_clean)
            if len(unique_alleles) == 1:
                allele = list(unique_alleles)[0]
                if allele == ref_allele:
                    for pos in range(start, stop):
                        if pos < len(sample_seqs[sample][contig]):
                            sample_seqs[sample][contig][pos] = reference_seqs[contig][pos]
                elif len(allele) == 1 and (stop - start) == 1:
                    sample_seqs[sample][contig][start] = allele.upper()
            elif (stop - start) == 1:
                sample_seqs[sample][contig][start] = (
                    get_iupac_code(called_alleles_clean) if use_iupac else reference_seqs[contig][start]
                )

    vcf.close()

    output_files = []
    for sample in samples:
        output_path = os.path.join(output_folder, f"{sample}_consensus.fasta")
        records = []
        for contig, seq_list in sample_seqs[sample].items():
            records.append(
                SeqRecord(
                    Seq("".join(seq_list)),
                    id=contig,
                    description=f"Consensus from VCF, sample={sample}, GQ>={min_gq}, DP>={min_dp}, IUPAC={use_iupac}",
                )
            )
        SeqIO.write(records, output_path, "fasta")
        output_files.append(output_path)
    return output_files


def main() -> None:
    reference_genome = (SCRIPT_DIR / REFERENCE_GENOME).resolve()
    vcf_dir = (SCRIPT_DIR / VCF_DIRECTORY).resolve()
    output_dir = (SCRIPT_DIR / OUTPUT_DIRECTORY).resolve()

    ensure_exists(reference_genome, "参考基因组文件")
    ensure_exists(vcf_dir, "VCF 输入目录")
    ensure_dir(output_dir)

    reference_seqs = parse_reference_fasta(str(reference_genome))
    vcf_patterns = [
        os.path.join(str(vcf_dir), "*.vcf.gz"),
        os.path.join(str(vcf_dir), "*.vcf"),
        os.path.join(str(vcf_dir), "*.g.vcf.gz"),
        os.path.join(str(vcf_dir), "*.g.vcf"),
    ]
    vcf_files = sorted(set(path for pattern in vcf_patterns for path in glob.glob(pattern)))
    if not vcf_files:
        raise FileNotFoundError(f"在 {vcf_dir} 中未找到 VCF 文件。")

    all_output_files: list[str] = []
    with ProcessPoolExecutor(max_workers=PARALLEL_SAMPLES) as executor:
        futures = {
            executor.submit(
                vcf_to_consensus,
                vcf_path,
                reference_seqs,
                str(output_dir),
                USE_IUPAC,
                MIN_GQ,
                MIN_DP,
            ): vcf_path
            for vcf_path in vcf_files
        }
        for future in as_completed(futures):
            all_output_files.extend(future.result())

    LOGGER.info(f"一致性序列生成完成，共输出 {len(sorted(all_output_files))} 个文件。")


if __name__ == "__main__":
    main()
