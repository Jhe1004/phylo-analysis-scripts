#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import os
import shutil
import subprocess

from Bio import SeqIO


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/1_mafft/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/12_WGD_detection/phypart/1_mafft/output"

CONDA_ENV_NAME = "trinity_env"
MAFFT_EXECUTABLE_NAME = "mafft"

PROCESS_COUNT = max(1, min(8, os.cpu_count() or 1))
THREADS_PER_PROCESS = 1
JUST_ALIGN_CDS_DIRECTLY = False
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def find_executable(executable_name: str) -> Path:
    local_candidate = SCRIPT_DIR / "dependencies" / "bin" / executable_name
    if local_candidate.exists():
        return local_candidate.resolve()
    conda_candidate = Path.home() / "miniconda3" / "envs" / CONDA_ENV_NAME / "bin" / executable_name
    if conda_candidate.exists():
        return conda_candidate.resolve()
    system_candidate = shutil.which(executable_name)
    if system_candidate:
        return Path(system_candidate).resolve()
    raise FileNotFoundError(f"未找到可执行文件 {executable_name}")


def aa_to_nt_alignment(aa_alignment_file: Path, cds_fasta_file: Path, output_file: Path) -> None:
    nt_dict = SeqIO.to_dict(SeqIO.parse(str(cds_fasta_file), "fasta"))
    aa_dict = SeqIO.to_dict(SeqIO.parse(str(aa_alignment_file), "fasta"))
    with open(output_file, "w", encoding="utf-8") as handle:
        for seq_id, aa_record in aa_dict.items():
            if seq_id not in nt_dict:
                raise KeyError(f"在 CDS 文件中找不到对应序列 ID: {seq_id}")
            nt_seq = str(nt_dict[seq_id].seq)
            if len(nt_seq) % 3 != 0:
                raise ValueError(f"CDS 长度不是 3 的倍数: {seq_id}")
            new_seq = []
            cursor = 0
            for aa_char in str(aa_record.seq):
                if aa_char == "-":
                    new_seq.append("---")
                elif aa_char == "*":
                    new_seq.append("---")
                    cursor += 3
                else:
                    new_seq.append(nt_seq[cursor : cursor + 3])
                    cursor += 3
            handle.write(f">{seq_id}\n{''.join(new_seq)}\n")


def run_command(cmd: list[str], output_file: Path) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"命令运行失败: {' '.join(cmd)}\n{result.stderr.strip()}")
    output_file.write_text(result.stdout, encoding="utf-8")


def process_family(task: tuple[Path, Path, Path, Path, bool, int]) -> str:
    pep_file, cds_file, pep_output, cds_output, just_align_cds, threads_per_process = task
    mafft_executable = find_executable(MAFFT_EXECUTABLE_NAME)
    run_command([str(mafft_executable), "--thread", str(threads_per_process), str(pep_file)], pep_output)
    if just_align_cds:
        run_command([str(mafft_executable), "--thread", str(threads_per_process), str(cds_file)], cds_output)
    else:
        aa_to_nt_alignment(pep_output, cds_file, cds_output)
    return pep_file.stem.replace("_pep", "")


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    pep_output_dir = output_dir / "pep_alignments"
    cds_output_dir = output_dir / "cds_alignments"
    ensure_directory(pep_output_dir)
    ensure_directory(cds_output_dir)

    pep_files = sorted(input_dir.glob("*_pep.fasta")) + sorted(input_dir.glob("*_pep.fas"))
    if not pep_files:
        raise FileNotFoundError(f"未在 {input_dir} 中找到 *_pep.fasta 文件")

    tasks = []
    for pep_file in pep_files:
        family_prefix = pep_file.name.rsplit("_pep", 1)[0]
        cds_file = input_dir / f"{family_prefix}_cds.fasta"
        if not cds_file.exists():
            cds_file = input_dir / f"{family_prefix}_cds.fas"
        if not cds_file.exists():
            raise FileNotFoundError(f"找不到与 {pep_file.name} 配套的 CDS 文件")
        tasks.append(
            (
                pep_file,
                cds_file,
                pep_output_dir / f"{family_prefix}_pep_maffted.fasta",
                cds_output_dir / f"{family_prefix}_cds_maffted.fasta",
                JUST_ALIGN_CDS_DIRECTLY,
                THREADS_PER_PROCESS,
            )
        )

    if PROCESS_COUNT > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=PROCESS_COUNT) as executor:
            results = list(executor.map(process_family, tasks))
    else:
        results = [process_family(task) for task in tasks]

    summary_file = output_dir / "mafft_summary.tsv"
    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write("family\tstatus\n")
        for family_name in results:
            handle.write(f"{family_name}\tsuccess\n")

    print(f"完成比对家族数: {len(results)}")
    print(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()
