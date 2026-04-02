#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path
import os
import shutil
import subprocess

from ete3 import Tree


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

CONDA_ENV_NAME = "trinity_env"
TREEPL_EXECUTABLE_NAME = "treePL"
TREEPL_LD_LIBRARY_PATH = "/usr/local/lib64"

CALIBRATED_AGE = 71.7
INGROUPS_FILE = "ingroupspecies.txt"
OUTGROUPS_FILE = "outgroupspecies.txt"
TREE_SUBDIRECTORY = "trees"
ALIGNMENT_SUBDIRECTORY = "alignments"

PROCESS_COUNT = max(1, min(12, os.cpu_count() or 1))
RUN_OUTPUT_PREFIX = "treepl_run_"

NTHREADS = 2
OPT = 3
OPTAD = 3
OPTCVAD = 5
CVITER = 5
CVSIMANITER = 10000
CVSTART = 10000
CVSTOP = 0.00000001
CVMULTSTEP = 0.1
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


def read_species_list(file_path: Path) -> list[str]:
    return [line.strip() for line in file_path.read_text(encoding="utf-8").splitlines() if line.strip()]


def get_species_in_clade(node) -> list[str]:
    result = []
    for leaf in node.get_leaf_names():
        species_name = leaf.split("++")[0]
        if species_name not in result:
            result.append(species_name)
    return result


def generate_calibration_list(tree_file: Path, ingroups_list: list[str], outgroups_list: list[str]) -> list[str]:
    tree = Tree(str(tree_file))

    def judge_outgroup_side(node):
        if len(node.children) != 2:
            return None
        left_species = get_species_in_clade(node.children[0])
        right_species = get_species_in_clade(node.children[1])
        left_has_outgroup = any(species in left_species for species in outgroups_list)
        right_has_outgroup = any(species in right_species for species in outgroups_list)
        if left_has_outgroup and right_has_outgroup:
            return None
        if left_has_outgroup and all(species in outgroups_list for species in left_species):
            return 0
        if right_has_outgroup and all(species in outgroups_list for species in right_species):
            return 1
        return None

    def is_large_ingroup_clade(node):
        species_in_node = get_species_in_clade(node)
        if any(species not in ingroups_list for species in species_in_node):
            return False
        ingroup_all = []
        for species in get_species_in_clade(tree):
            if species in ingroups_list and species not in ingroup_all:
                ingroup_all.append(species)
        if not ingroup_all:
            return False
        return len(species_in_node) / len(ingroup_all) > (2 / 3)

    calibrations = []
    counter = 0
    for node in tree.traverse():
        if not node.children:
            continue
        outgroup_side = judge_outgroup_side(node)
        if outgroup_side is None:
            continue
        candidate = node.children[1] if outgroup_side == 0 else node.children[0]
        if not is_large_ingroup_clade(candidate):
            continue
        calibrations.append(f"mrca = a{counter} {' '.join(node.get_leaf_names())}\n")
        calibrations.append(f"min = a{counter} {CALIBRATED_AGE}\n")
        calibrations.append(f"max = a{counter} {CALIBRATED_AGE + 0.01}\n")
        counter += 1
    return calibrations


def find_alignment_file(tree_file: Path, alignment_dir: Path) -> Path:
    alignment_name = tree_file.name
    if alignment_name.startswith("RAxML_bestTree."):
        alignment_name = alignment_name.replace("RAxML_bestTree.", "", 1)
    if alignment_name.startswith("RAxML_bipartitions."):
        alignment_name = alignment_name.replace("RAxML_bipartitions.", "", 1)
    if alignment_name.endswith("_r.tree"):
        alignment_name = alignment_name[:-7] + ".fas"
    else:
        alignment_name = alignment_name[:-5] + ".fas"
    alignment_path = alignment_dir / alignment_name
    if not alignment_path.exists():
        raise FileNotFoundError(f"找不到对应比对文件: {alignment_path}")
    return alignment_path


def read_sequence_length(alignment_file: Path) -> int:
    with open(alignment_file, "r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if text and not text.startswith(">"):
                return len(text)
    raise ValueError(f"未在 {alignment_file} 中找到有效序列")


def write_treepl_conf(
    tree_file: Path,
    output_conf: Path,
    alignment_file: Path,
    calibrations: list[str],
    smooth_value: float,
    cv_mode: bool,
) -> Path:
    cvoutfile = Path(str(output_conf) + ".cvout")
    with open(output_conf, "w", encoding="utf-8") as handle:
        handle.write("[Input files containing the ML trees]\n")
        handle.write(f"treefile = {tree_file.resolve()}\n")
        handle.write("[General commands]\n")
        handle.write(f"nthreads = {NTHREADS}\n")
        handle.write(f"smooth = {smooth_value}\n")
        handle.write("thorough\n")
        handle.write("log_pen\n")
        handle.write(f"numsites = {read_sequence_length(alignment_file)}\n")
        handle.write("[Calibrations]\n")
        for line in calibrations:
            handle.write(line)
        handle.write("[Best optimisation parameters]\n")
        handle.write(f"opt = {OPT}\n")
        handle.write("moredetail\n")
        handle.write(f"optad = {OPTAD}\n")
        handle.write("moredetailad\n")
        handle.write(f"optcvad = {OPTCVAD}\n")
        if cv_mode:
            handle.write("cv\n")
            handle.write("randomcv\n")
            handle.write(f"cviter = {CVITER}\n")
            handle.write(f"cvsimaniter = {CVSIMANITER}\n")
            handle.write(f"cvstart = {CVSTART}\n")
            handle.write(f"cvstop = {CVSTOP}\n")
            handle.write(f"cvmultstep = {CVMULTSTEP}\n")
            handle.write(f"cvoutfile = {cvoutfile}\n")
        handle.write(f"outfile = {output_conf}.newtree\n")
    return cvoutfile


def run_treepl(treepl_executable: Path, conf_file: Path, log_file: Path) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    if TREEPL_LD_LIBRARY_PATH:
        current = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = f"{TREEPL_LD_LIBRARY_PATH}:{current}" if current else TREEPL_LD_LIBRARY_PATH
    result = subprocess.run(
        [str(treepl_executable), str(conf_file)],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    log_file.write_text((result.stdout or "") + (result.stderr or ""), encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"treePL 运行失败: {conf_file.name}\n{result.stderr.strip()}")
    return result


def parse_cv_results(cvoutfile: Path) -> float:
    results = []
    with open(cvoutfile, "r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if text.startswith("chisq:"):
                left, right = text.split(")", maxsplit=1)
                smooth_text = left.split("(", maxsplit=1)[1].strip()
                chisq_text = right.strip()
                results.append((float(smooth_text), float(chisq_text)))
    if not results:
        raise ValueError(f"未在 {cvoutfile} 中解析到 cross-validation 结果")
    best_smooth, _ = min(results, key=lambda item: item[1])
    return best_smooth


def process_tree(task: tuple[Path, Path, Path, list[str], list[str], Path]) -> tuple[str, str, str]:
    tree_file, alignment_dir, run_dir, ingroups_list, outgroups_list, treepl_executable = task
    try:
        calibrations = generate_calibration_list(tree_file, ingroups_list, outgroups_list)
        if not calibrations:
            return tree_file.name, "skip", "no_calibration"
        alignment_file = find_alignment_file(tree_file, alignment_dir)
        log_dir = run_dir / "logs"
        ensure_directory(log_dir)
        output_conf = run_dir / f"{tree_file.stem}.conf"
        cv_conf = run_dir / f"{output_conf.name}.cv.conf"
        cvoutfile = write_treepl_conf(tree_file, cv_conf, alignment_file, calibrations, CVSTART, True)
        run_treepl(treepl_executable, cv_conf, log_dir / f"{cv_conf.name}.log")
        best_smooth = parse_cv_results(cvoutfile)
        write_treepl_conf(tree_file, output_conf, alignment_file, calibrations, best_smooth, False)
        run_treepl(treepl_executable, output_conf, log_dir / f"{output_conf.name}.log")
        return tree_file.name, "success", str(best_smooth)
    except Exception as exc:
        return tree_file.name, "error", str(exc).replace("\n", " ")


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    ensure_directory(output_dir)

    tree_dir = input_dir / TREE_SUBDIRECTORY
    alignment_dir = input_dir / ALIGNMENT_SUBDIRECTORY
    ingroups_list = read_species_list(input_dir / INGROUPS_FILE)
    outgroups_list = read_species_list(input_dir / OUTGROUPS_FILE)
    treepl_executable = find_executable(TREEPL_EXECUTABLE_NAME)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = output_dir / f"{RUN_OUTPUT_PREFIX}{timestamp}"
    ensure_directory(run_dir)

    tree_files = sorted(tree_dir.glob("*.tree"))
    if not tree_files:
        raise FileNotFoundError(f"未在 {tree_dir} 中找到 .tree 文件")

    tasks = [(tree_file, alignment_dir, run_dir, ingroups_list, outgroups_list, treepl_executable) for tree_file in tree_files]
    if PROCESS_COUNT > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=PROCESS_COUNT) as executor:
            results = list(executor.map(process_tree, tasks))
    else:
        results = [process_tree(task) for task in tasks]

    summary_file = run_dir / "best_smooth_summary.tsv"
    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write("tree\tstatus\tbest_smooth_or_message\n")
        for tree_name, status, message in results:
            handle.write(f"{tree_name}\t{status}\t{message}\n")

    print(f"运行目录: {run_dir}")
    print(f"处理树数: {len(results)}")


if __name__ == "__main__":
    main()
