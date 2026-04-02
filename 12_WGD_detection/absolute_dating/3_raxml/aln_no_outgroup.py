#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import os
import shutil
import subprocess


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

CONDA_ENV_NAME = "trinity_env"
RAXML_EXECUTABLE_NAME = "raxmlHPC-PTHREADS"

INPUT_EXTENSION = ".fas"
PROCESS_COUNT = max(1, min(12, os.cpu_count() or 1))
THREADS_PER_PROCESS = 1
MODEL = "GTRGAMMA"
RANDOM_SEED = 12345
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


def run_single_alignment(
    raxml_executable: Path,
    alignment_file: Path,
    output_dir: Path,
    log_dir: Path,
) -> tuple[str, int]:
    name = alignment_file.stem
    cmd = [
        str(raxml_executable),
        "-T",
        str(THREADS_PER_PROCESS),
        "-n",
        name,
        "-s",
        str(alignment_file),
        "-m",
        MODEL,
        "-p",
        str(RANDOM_SEED),
        "-w",
        str(output_dir.resolve()),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (log_dir / f"{name}.log").write_text((result.stdout or "") + (result.stderr or ""), encoding="utf-8")
    return alignment_file.name, result.returncode


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    log_dir = output_dir / "logs"
    ensure_directory(output_dir)
    ensure_directory(log_dir)

    alignment_files = sorted(input_dir.glob(f"*{INPUT_EXTENSION}"))
    if not alignment_files:
        raise FileNotFoundError(f"未在 {input_dir} 中找到 {INPUT_EXTENSION} 文件")
    raxml_executable = find_executable(RAXML_EXECUTABLE_NAME)
    print(f"使用 RAxML: {raxml_executable}")

    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=PROCESS_COUNT) as executor:
        futures = {
            executor.submit(run_single_alignment, raxml_executable, file_path, output_dir, log_dir): file_path
            for file_path in alignment_files
        }
        for future in as_completed(futures):
            file_name, return_code = future.result()
            if return_code != 0:
                failures.append(file_name)

    summary_file = output_dir / "raxml_summary.tsv"
    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write("file\tstatus\n")
        for alignment_file in alignment_files:
            status = "failed" if alignment_file.name in failures else "success"
            handle.write(f"{alignment_file.name}\t{status}\n")

    if failures:
        raise RuntimeError(f"以下比对建树失败: {', '.join(failures)}")
    print(f"完成建树数: {len(alignment_files)}")


if __name__ == "__main__":
    main()
