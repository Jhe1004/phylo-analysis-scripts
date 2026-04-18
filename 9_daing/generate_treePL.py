#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
单树 single-file treePL dating 流程。
"""

from __future__ import annotations

from pathlib import Path
import os
import shutil
import subprocess

from ete3 import Tree


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/9_daing/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/9_daing/output"

CONDA_ENV_NAME = "trinity_env"
TREEPL_EXECUTABLE_NAME = "treePL"
TREEPL_LD_LIBRARY_PATH = "/usr/local/lib64"

TREE_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/9_daing/input/coalescent_ml.treefile"
ALIGNMENT_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/9_daing/input/result0.1.fas"

ROOT_AGE = None
CALIBRATION_TEXT = """
# 每个节点占 3 行，可连续写多个节点。
# 第 1 行：节点名 + 一个子树 Newick 片段
# 第 2 行：节点名max + 最大年龄
# 第 3 行：节点名min + 最小年龄
#
# 可直接参考下面这个示例，按同样格式继续往下写：
# node1 ((Species_A:0.1,Species_B:0.2):0.3,Species_C:0.4);
# node1max 72
# node1min 68
#
# node2 (Species_D:0.1,(Species_E:0.2,Species_F:0.3):0.4);
# node2max 55
# node2min 50
"""

USE_RANDOMCV = True
CV_ITER = 5
CV_SIMAN_ITER = 100000
CV_START = 10000
CV_STOP = 0.00001
CV_MULTSTEP = 0.1

NTHREADS = 12
THOROUGH = True
LOG_PEN = True
OPT = 3
OPTAD = 3
OPTCVAD = 5
MOREDETAIL = True
MOREDETAILAD = True
# ==============================


OUTPUT_CONF_FILE = "dating.conf"


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

    raise FileNotFoundError(
        f"未找到可执行文件 {executable_name}。查找顺序为 dependencies/bin -> conda 环境 -> PATH。"
    )


def parse_calibration_text(calibration_text: str) -> list[str]:
    lines = []
    for raw_line in calibration_text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        lines.append(line)
    return lines


def subtree_to_leaf_names(subtree_text: str) -> list[str]:
    subtree_text = subtree_text.strip()
    if not subtree_text:
        raise ValueError("校准子树不能为空")
    if not subtree_text.endswith(";"):
        subtree_text += ";"
    subtree = Tree(subtree_text, format=1)
    leaves = subtree.get_leaf_names()
    if len(leaves) < 2:
        raise ValueError(f"校准子树至少需要 2 个叶节点: {subtree_text}")
    return leaves


def build_text_calibrations(tree_file: Path, calibration_text: str) -> list[str]:
    tree = Tree(str(tree_file), format=1)
    raw_lines = parse_calibration_text(calibration_text)
    if not raw_lines:
        return []
    if len(raw_lines) % 3 != 0:
        raise ValueError(
            "CALIBRATION_TEXT 的有效行数必须是 3 的倍数：每个节点需要 node / nodemax / nodemin 三行"
        )

    calibration_lines: list[str] = []
    used_names: set[str] = set()
    for index in range(0, len(raw_lines), 3):
        node_line = raw_lines[index]
        bound_line_1 = raw_lines[index + 1]
        bound_line_2 = raw_lines[index + 2]

        node_parts = node_line.split(maxsplit=1)
        if len(node_parts) != 2:
            raise ValueError(f"节点定义行格式错误，应为 'node_name <subtree_newick>': {node_line}")
        node_name, subtree_text = node_parts
        if node_name in used_names:
            raise ValueError(f"重复的节点名: {node_name}")
        used_names.add(node_name)

        bounds = {}
        for bound_line in (bound_line_1, bound_line_2):
            bound_parts = bound_line.split(maxsplit=1)
            if len(bound_parts) != 2:
                raise ValueError(f"年龄定义行格式错误: {bound_line}")
            key, value = bound_parts
            if key == f"{node_name}max":
                bounds["max"] = float(value)
            elif key == f"{node_name}min":
                bounds["min"] = float(value)
            else:
                raise ValueError(f"年龄定义行与节点名不匹配: {bound_line}")

        if "max" not in bounds or "min" not in bounds:
            raise ValueError(f"节点 {node_name} 缺少 max 或 min 年龄")
        if bounds["min"] > bounds["max"]:
            raise ValueError(f"节点 {node_name} 的最小年龄大于最大年龄")

        selected_leaves = subtree_to_leaf_names(subtree_text)
        try:
            mrca_node = tree.get_common_ancestor(selected_leaves)
        except Exception as exc:
            raise ValueError(
                f"节点 {node_name} 的子树在目标树中找不到，请检查物种名是否完全一致"
            ) from exc

        mrca_leaves = mrca_node.get_leaf_names()
        if set(mrca_leaves) != set(selected_leaves):
            raise ValueError(
                f"节点 {node_name} 对应的子树在目标树里不是单系。"
                f" 给定叶节点数={len(selected_leaves)}，实际 MRCA 覆盖数={len(mrca_leaves)}"
            )

        calibration_lines.append(f"mrca = {node_name} {' '.join(selected_leaves)}\n")
        calibration_lines.append(f"max = {node_name} {bounds['max']}\n")
        calibration_lines.append(f"min = {node_name} {bounds['min']}\n")

    return calibration_lines


def generate_calibration_list(tree_file: Path, calibration_text: str, root_age: float | None) -> list[str]:
    text_calibrations = build_text_calibrations(tree_file, calibration_text)
    if text_calibrations:
        return text_calibrations

    if root_age is None:
        raise ValueError("未检测到 CALIBRATION_TEXT 中的校准节点，同时也没有提供根节点年龄")

    tree = Tree(str(tree_file), format=1)
    all_leaves = tree.get_leaf_names()
    if not all_leaves:
        raise ValueError(f"树文件中没有叶节点: {tree_file}")

    return [
        f"mrca = root_node {' '.join(all_leaves)}\n",
        f"fixage = root_node {root_age}\n",
    ]


def read_sequence_length(alignment_file: Path) -> int:
    with open(alignment_file, "r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith(">"):
                continue
            return len(text)
    raise ValueError(f"未在 {alignment_file} 中找到有效序列行")


def write_general_settings(handle, seq_len: int, smooth_value: float) -> None:
    handle.write("[General commands]\n")
    handle.write(f"nthreads = {NTHREADS}\n")
    handle.write(f"smooth = {smooth_value}\n")
    if THOROUGH:
        handle.write("thorough\n")
    if LOG_PEN:
        handle.write("log_pen\n")
    handle.write(f"numsites = {seq_len}\n")


def write_optimization_settings(handle) -> None:
    handle.write("[Best optimisation parameters]\n")
    handle.write(f"opt = {OPT}\n")
    if MOREDETAIL:
        handle.write("moredetail\n")
    handle.write(f"optad = {OPTAD}\n")
    if MOREDETAILAD:
        handle.write("moredetailad\n")
    handle.write(f"optcvad = {OPTCVAD}\n")


def write_cv_settings(handle, cvoutfile: Path) -> None:
    handle.write("cv\n")
    if USE_RANDOMCV:
        handle.write("randomcv\n")
    handle.write(f"cviter = {CV_ITER}\n")
    handle.write(f"cvsimaniter = {CV_SIMAN_ITER}\n")
    handle.write(f"cvstart = {CV_START}\n")
    handle.write(f"cvstop = {CV_STOP}\n")
    handle.write(f"cvmultstep = {CV_MULTSTEP}\n")
    handle.write(f"cvoutfile = {cvoutfile}\n")


def generate_treepl_conf(
    tree_file: Path,
    alignment_file: Path,
    calibration_text: str,
    root_age: float | None,
    output_conf: Path,
    smooth_value: float,
    cv_mode: bool = False,
) -> Path:
    seq_len = read_sequence_length(alignment_file)
    calib_lines = generate_calibration_list(tree_file, calibration_text, root_age)
    cvoutfile = output_conf.with_suffix(output_conf.suffix + ".cvout")

    with open(output_conf, "w", encoding="utf-8") as handle:
        handle.write("[Input files containing the ML trees]\n")
        handle.write(f"treefile = {tree_file}\n")
        write_general_settings(handle, seq_len, smooth_value)
        handle.write("[Calibrations]\n")
        for line in calib_lines:
            handle.write(line)
        write_optimization_settings(handle)
        if cv_mode:
            write_cv_settings(handle, cvoutfile)
        handle.write(f"outfile = {output_conf}.newtree\n")

    return cvoutfile


def parse_cv_results(cvoutfile: Path) -> tuple[float, list[tuple[float, float]]]:
    results = []
    with open(cvoutfile, "r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text.startswith("chisq:"):
                continue
            left, right = text.split(")", maxsplit=1)
            smooth_text = left.split("(", maxsplit=1)[1].strip()
            chisq_text = right.strip()
            results.append((float(smooth_text), float(chisq_text)))
    if not results:
        raise ValueError(f"未在 {cvoutfile} 中解析到 cross-validation 结果")
    best_smooth, _ = min(results, key=lambda item: item[1])
    return best_smooth, results


def run_treepl(treepl_bin: Path, conf_file: Path, log_file: Path) -> subprocess.CompletedProcess:
    cmd = [str(treepl_bin), str(conf_file)]
    env = dict(os.environ)
    if TREEPL_LD_LIBRARY_PATH:
        current_ld_library_path = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = (
            f"{TREEPL_LD_LIBRARY_PATH}:{current_ld_library_path}"
            if current_ld_library_path
            else TREEPL_LD_LIBRARY_PATH
        )
    try:
        completed = subprocess.run(
            cmd,
            check=True,
            text=True,
            capture_output=True,
            env=env,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"找不到 treePL 可执行文件: {treepl_bin}") from exc
    except subprocess.CalledProcessError as exc:
        message = exc.stderr.strip() or exc.stdout.strip() or str(exc)
        raise RuntimeError(f"treePL 运行失败: {' '.join(cmd)}\n{message}") from exc

    with open(log_file, "w", encoding="utf-8") as handle:
        if completed.stdout:
            handle.write(completed.stdout)
            if not completed.stdout.endswith("\n"):
                handle.write("\n")
        if completed.stderr:
            handle.write(completed.stderr)
            if not completed.stderr.endswith("\n"):
                handle.write("\n")
    return completed


def run_cv_and_dating(
    tree_file: Path,
    alignment_file: Path,
    calibration_text: str,
    root_age: float | None,
    output_conf: Path,
    treepl_bin: Path,
    log_dir: Path,
) -> tuple[Path, Path, float, list[tuple[float, float]]]:
    cv_conf = output_conf.with_suffix(output_conf.suffix + ".cv.conf")
    cvoutfile = generate_treepl_conf(
        tree_file=tree_file,
        alignment_file=alignment_file,
        calibration_text=calibration_text,
        root_age=root_age,
        output_conf=cv_conf,
        smooth_value=CV_START,
        cv_mode=True,
    )
    cv_run = run_treepl(treepl_bin, cv_conf, log_dir / f"{cv_conf.name}.log")
    if not cvoutfile.exists():
        details = cv_run.stderr.strip() or cv_run.stdout.strip() or "treePL 未生成 cvoutfile"
        raise RuntimeError(
            "treePL 已执行 cv 配置，但没有生成 cvoutfile。"
            f"\n配置文件: {cv_conf}"
            f"\n期望的 cvoutfile: {cvoutfile}"
            f"\n程序输出:\n{details}"
        )

    best_smooth, cv_results = parse_cv_results(cvoutfile)

    final_cvoutfile = generate_treepl_conf(
        tree_file=tree_file,
        alignment_file=alignment_file,
        calibration_text=calibration_text,
        root_age=root_age,
        output_conf=output_conf,
        smooth_value=best_smooth,
        cv_mode=False,
    )
    if final_cvoutfile.exists():
        final_cvoutfile.unlink()
    run_treepl(treepl_bin, output_conf, log_dir / f"{output_conf.name}.log")
    return cv_conf, cvoutfile, best_smooth, cv_results


def main() -> None:
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    log_dir = output_dir / "logs"
    ensure_directory(output_dir)
    ensure_directory(log_dir)

    treepl_bin = find_executable(TREEPL_EXECUTABLE_NAME)
    tree_in = input_dir / TREE_FILE
    alignment_in = input_dir / ALIGNMENT_FILE
    output_conf = output_dir / OUTPUT_CONF_FILE

    print(f"使用 treePL: {treepl_bin}")

    if not tree_in.exists():
        raise FileNotFoundError(f"找不到树文件: {tree_in}")
    if not alignment_in.exists():
        raise FileNotFoundError(f"找不到比对文件: {alignment_in}")

    calibration_lines = parse_calibration_text(CALIBRATION_TEXT)
    root_age = ROOT_AGE
    if not calibration_lines:
        if root_age is None:
            raise ValueError("未设置 CALIBRATION_TEXT 时，必须在脚本顶部设置 ROOT_AGE")
        try:
            root_age = float(root_age)
        except ValueError as exc:
            raise ValueError(f"ROOT_AGE 必须是数值: {root_age}") from exc
        if root_age <= 0:
            raise ValueError(f"根节点校准年龄必须大于 0: {root_age}")

    cv_conf, cvoutfile, best_smooth, cv_results = run_cv_and_dating(
        tree_file=tree_in,
        alignment_file=alignment_in,
        calibration_text=CALIBRATION_TEXT,
        root_age=root_age,
        output_conf=output_conf,
        treepl_bin=treepl_bin,
        log_dir=log_dir,
    )

    summary_file = output_dir / "dating_summary.txt"
    with open(summary_file, "w", encoding="utf-8") as handle:
        handle.write(f"treePL = {treepl_bin}\n")
        handle.write(f"tree_file = {tree_in}\n")
        handle.write(f"alignment_file = {alignment_in}\n")
        handle.write(f"cv_conf = {cv_conf}\n")
        handle.write(f"cvoutfile = {cvoutfile}\n")
        handle.write(f"best_smooth = {best_smooth}\n")
        handle.write(f"final_conf = {output_conf}\n")
        handle.write(f"final_tree = {output_conf}.newtree\n")
        handle.write("cv_results =\n")
        for smooth_value, chisq_value in cv_results:
            handle.write(f"  smooth={smooth_value}\tchisq={chisq_value}\n")

    print(f"CV 配置文件: {cv_conf}")
    print(f"CV 结果文件: {cvoutfile}")
    print(f"最佳 smooth: {best_smooth}")
    print(f"最终配置文件: {output_conf}")
    print(f"结果摘要: {summary_file}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"错误: {exc}")
        raise SystemExit(1)
