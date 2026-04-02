#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量运行 treePL，并支持多任务并行。
"""

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
TREEPL_EXECUTABLE_NAME = "treePL"

CONFIG_SUBDIRECTORY = "configs"
PROCESS_COUNT = max(1, min(6, os.cpu_count() or 1))
SKIP_COMPLETED = True
# ==============================


def resolve_path(relative_path: str) -> Path:
    return (SCRIPT_DIR / relative_path).resolve()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def find_executable(executable_name: str) -> Path:
    local_candidate = SCRIPT_DIR / "dependencies" / "bin" / executable_name
    if local_candidate.exists():
        return local_candidate.resolve()

    conda_prefix = Path.home() / "miniconda3" / "envs" / CONDA_ENV_NAME / "bin" / executable_name
    if conda_prefix.exists():
        return conda_prefix.resolve()

    system_candidate = shutil.which(executable_name)
    if system_candidate:
        return Path(system_candidate).resolve()

    raise FileNotFoundError(
        f"未找到可执行文件 {executable_name}。查找顺序为 dependencies/bin -> conda 环境 -> PATH。"
    )


def copy_input_workspace(input_dir: Path, output_workspace: Path) -> None:
    if output_workspace.exists():
        shutil.rmtree(output_workspace)
    shutil.copytree(input_dir, output_workspace)


def parse_outfile_from_conf(conf_file: Path) -> Path | None:
    with open(conf_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip().startswith("outfile"):
                _, value = line.split("=", 1)
                return (conf_file.parent / value.strip()).resolve()
    return None


def run_single_treepl(
    conf_file: Path,
    treepl_executable: Path,
    results_dir: Path,
    log_dir: Path,
) -> tuple[str, int]:
    output_tree = parse_outfile_from_conf(conf_file)
    copied_output = results_dir / f"{conf_file.stem}.newtree"
    if SKIP_COMPLETED and copied_output.exists():
        return conf_file.name, 0

    log_file = log_dir / f"{conf_file.stem}.log"
    with open(log_file, "w", encoding="utf-8") as handle:
        process = subprocess.run(
            [str(treepl_executable), conf_file.name],
            cwd=str(conf_file.parent),
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )

    if process.returncode == 0 and output_tree and output_tree.exists():
        shutil.copy2(output_tree, copied_output)

    return conf_file.name, process.returncode


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    workspace_dir = output_dir / "workspace"
    results_dir = output_dir / "dated_trees"
    log_dir = output_dir / "logs"

    ensure_directory(output_dir)
    ensure_directory(results_dir)
    ensure_directory(log_dir)

    treepl_executable = find_executable(TREEPL_EXECUTABLE_NAME)
    print(f"使用 treePL: {treepl_executable}")

    copy_input_workspace(input_dir, workspace_dir)
    conf_dir = workspace_dir / CONFIG_SUBDIRECTORY
    conf_files = sorted(conf_dir.glob("*.conf"))
    if not conf_files:
        raise FileNotFoundError(f"未找到配置文件: {conf_dir}")

    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=PROCESS_COUNT) as executor:
        futures = {
            executor.submit(run_single_treepl, conf_file, treepl_executable, results_dir, log_dir): conf_file
            for conf_file in conf_files
        }
        for future in as_completed(futures):
            conf_name, return_code = future.result()
            if return_code != 0:
                failures.append(conf_name)
                print(f"失败: {conf_name}")
            else:
                print(f"完成: {conf_name}")

    summary_file = output_dir / "treepl_run_summary.tsv"
    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write("conf_file\tstatus\n")
        for conf_file in conf_files:
            status = "failed" if conf_file.name in failures else "success"
            handle.write(f"{conf_file.name}\t{status}\n")

    if failures:
        raise RuntimeError(f"以下配置运行失败: {', '.join(failures)}")

    print(f"全部完成。结果目录: {output_dir}")


if __name__ == "__main__":
    main()
