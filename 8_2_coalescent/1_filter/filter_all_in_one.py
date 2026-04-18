from __future__ import annotations

import glob
import math
import os
import shutil
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import dendropy
from Bio import SeqIO


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_2_coalescent/1_filter/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_2_coalescent/1_filter/output"
CONDA_ENV_NAME = "trinity_env"

FASTA_GLOB = "ortho*_cds_maffted.fas"
TREE_GLOB = "RAxML_bipartitions.ortho*_cds_maffted"
FASTA_NAME_TEMPLATE = "{basename}.fas"
TREE_NAME_TEMPLATE = "RAxML_bipartitions.{basename}"

AUTO_PREPARE_ASTRAL_INPUTS = True
ASTRAL_INPUT_RELATIVE_DIR = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_2_coalescent/2_0_astral/input"
ASTRAL_TREE_FILE_GLOB = "RAxML_bipartitions.*"

ENABLE_PARALLEL = True
MAX_WORKERS = min(32, max(1, (os.cpu_count() or 1) - 1))

RUN_MISSINGDATA_FILTER = True
RUN_BOOTSTRAP_FILTER = True
RUN_SEQ_LENGTH_FILTER = True
RUN_PI_SITE_FILTER = True

MISSING_DATA_THRESHOLDS = [0.1, 0.2]
BOOTSTRAP_THRESHOLDS = [60, 70]
SEQ_LENGTH_THRESHOLDS = [1000, 1500]
PI_SITE_THRESHOLDS = [500, 1000]

CLEAR_OUTPUT_DIRS_BEFORE_RUN = True
COPY_FASTA_EVEN_IF_TREE_MISSING = True
VALID_BASES = {"A", "T", "C", "G"}
# ============================================================


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = Path(INPUT_DIRECTORY).resolve()
OUTPUT_ROOT = Path(OUTPUT_DIRECTORY).resolve()
ASTRAL_INPUT_DIR = Path(ASTRAL_INPUT_RELATIVE_DIR).resolve()


@dataclass
class LocusMetrics:
    basename: str
    fasta_path: Path
    tree_path: Path | None
    taxon_count: int
    seq_length: int
    missing_data_ratio: float
    avg_bootstrap: float | None
    pi_sites: int


@dataclass
class MetricTasks:
    need_bootstrap: bool
    need_seq_length: bool
    need_pi_sites: bool


def ensure_dir(path: Path) -> None:
    if CLEAR_OUTPUT_DIRS_BEFORE_RUN and path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def discover_fasta_files() -> list[Path]:
    fasta_paths = sorted(Path(p) for p in glob.glob(str(INPUT_DIR / FASTA_GLOB)))
    return [p for p in fasta_paths if p.is_file()]


def infer_tree_path(basename: str) -> Path:
    return INPUT_DIR / TREE_NAME_TEMPLATE.format(basename=basename)


def fasta_basename_from_path(fasta_path: Path) -> str:
    suffix = fasta_path.suffix
    return fasta_path.name[: -len(suffix)] if suffix else fasta_path.name


def compute_average_bootstrap(tree_path: Path) -> float | None:
    try:
        tree = dendropy.Tree.get(path=str(tree_path), schema="newick")
    except Exception as exc:
        print(f"[WARN] 树文件读取失败，跳过 bootstrap 统计: {tree_path.name} ({exc})")
        return None
    values = []
    for node in tree.preorder_node_iter():
        if node.label is None:
            continue
        try:
            values.append(float(node.label))
        except ValueError:
            continue
    return sum(values) / len(values) if values else None


def compute_pi_sites(records, source_name: str) -> tuple[int, int]:
    if not records:
        return 0, 0
    seq_lengths = {len(record.seq) for record in records}
    seq_length = min(seq_lengths) if len(seq_lengths) != 1 else max(seq_lengths)
    if len(seq_lengths) != 1:
        print(f"[WARN] 序列长度不一致，按最短公共范围统计 PI 位点: {source_name}")
    pi_sites = 0
    for i in range(seq_length):
        counts = Counter()
        for rec in records:
            base = str(rec.seq[i]).upper()
            if base in VALID_BASES:
                counts[base] += 1
        if len(counts) >= 2 and sum(1 for count in counts.values() if count >= 2) >= 2:
            pi_sites += 1
    return pi_sites, seq_length


def analyze_locus(task_input: tuple[Path, MetricTasks]) -> LocusMetrics | None:
    fasta_path, tasks = task_input
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        print(f"[WARN] 空 fasta 文件，跳过: {fasta_path.name}")
        return None
    basename = fasta_basename_from_path(fasta_path)
    tree_path = infer_tree_path(basename)
    tree_path_obj = tree_path if tree_path.exists() else None
    taxon_count = len(records)
    if tasks.need_seq_length or tasks.need_pi_sites:
        pi_sites, seq_length = compute_pi_sites(records, fasta_path.name)
        if not tasks.need_pi_sites:
            pi_sites = 0
    else:
        seq_length = len(records[0].seq)
        pi_sites = 0
    avg_bootstrap = compute_average_bootstrap(tree_path_obj) if tasks.need_bootstrap and tree_path_obj else None
    return LocusMetrics(
        basename=basename,
        fasta_path=fasta_path,
        tree_path=tree_path_obj,
        taxon_count=taxon_count,
        seq_length=seq_length if tasks.need_seq_length or tasks.need_pi_sites else 0,
        missing_data_ratio=0.0,
        avg_bootstrap=avg_bootstrap,
        pi_sites=pi_sites,
    )


def collect_metrics() -> list[LocusMetrics]:
    fasta_files = discover_fasta_files()
    if not fasta_files:
        raise FileNotFoundError(f"没有找到匹配的 fasta 文件。INPUT_DIR={INPUT_DIR}, FASTA_GLOB={FASTA_GLOB}")
    tasks = MetricTasks(
        need_bootstrap=RUN_BOOTSTRAP_FILTER,
        need_seq_length=RUN_SEQ_LENGTH_FILTER,
        need_pi_sites=RUN_PI_SITE_FILTER,
    )
    task_inputs = [(fasta_path, tasks) for fasta_path in fasta_files]
    if ENABLE_PARALLEL and len(task_inputs) > 1:
        worker_count = min(MAX_WORKERS, len(task_inputs))
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            metrics_candidates = list(executor.map(analyze_locus, task_inputs))
    else:
        metrics_candidates = [analyze_locus(task_input) for task_input in task_inputs]
    metrics_list = [metrics for metrics in metrics_candidates if metrics is not None]
    max_taxon_count = max(metrics.taxon_count for metrics in metrics_list)
    for metrics in metrics_list:
        metrics.missing_data_ratio = 1 - (metrics.taxon_count / max_taxon_count)
        if not (RUN_SEQ_LENGTH_FILTER or RUN_PI_SITE_FILTER):
            records = list(SeqIO.parse(str(metrics.fasta_path), "fasta"))
            metrics.seq_length = len(records[0].seq) if records else 0
    return metrics_list


def safe_threshold_label(value: float | int) -> str:
    if isinstance(value, float) and not value.is_integer():
        text = f"{value:.6f}".rstrip("0").rstrip(".")
    else:
        text = str(int(value)) if float(value).is_integer() else str(value)
    return text.replace(".", "p")


def copy_locus_files(locus: LocusMetrics, output_dir: Path) -> None:
    shutil.copyfile(locus.fasta_path, output_dir / locus.fasta_path.name)
    if locus.tree_path is not None and locus.tree_path.exists():
        shutil.copyfile(locus.tree_path, output_dir / locus.tree_path.name)
    elif not COPY_FASTA_EVEN_IF_TREE_MISSING:
        target = output_dir / locus.fasta_path.name
        if target.exists():
            target.unlink()


def export_by_thresholds(metrics_list, thresholds, metric_name, selector, value_formatter):
    metric_root = OUTPUT_ROOT / metric_name
    ensure_dir(metric_root)
    for threshold in thresholds:
        out_dir = metric_root / f"{metric_name}_{safe_threshold_label(threshold)}"
        ensure_dir(out_dir)
        passed = [locus for locus in metrics_list if selector(locus, threshold)]
        passed_basenames = {locus.basename for locus in passed}
        for locus in passed:
            copy_locus_files(locus, out_dir)
        summary_path = out_dir / "summary.tsv"
        with summary_path.open("w", encoding="utf-8") as handle:
            handle.write("basename\tfasta\ttree\ttaxon_count\tseq_length\tmissing_data_ratio\tavg_bootstrap\tpi_sites\tthreshold\tresult\n")
            for locus in metrics_list:
                result = "PASS" if locus.basename in passed_basenames else "FAIL"
                tree_name = locus.tree_path.name if locus.tree_path is not None else "NA"
                bootstrap_text = "NA" if locus.avg_bootstrap is None or math.isnan(locus.avg_bootstrap) else f"{locus.avg_bootstrap:.4f}"
                handle.write(
                    f"{locus.basename}\t{locus.fasta_path.name}\t{tree_name}\t{locus.taxon_count}\t{locus.seq_length}\t"
                    f"{locus.missing_data_ratio:.6f}\t{bootstrap_text}\t{locus.pi_sites}\t{value_formatter(threshold)}\t{result}\n"
                )


def normalize_tree_text(tree_text: str, source_name: str) -> str | None:
    normalized = tree_text.strip()
    if not normalized:
        return None
    if not normalized.endswith(";"):
        normalized += ";"
    return normalized


def prepare_astral_inputs() -> None:
    if not AUTO_PREPARE_ASTRAL_INPUTS:
        return
    ASTRAL_INPUT_DIR.mkdir(parents=True, exist_ok=True)
    for tree_file in glob.glob(str(INPUT_DIR / ASTRAL_TREE_FILE_GLOB)):
        source = Path(tree_file)
        target = ASTRAL_INPUT_DIR / source.name
        shutil.copyfile(source, target)


def main():
    if not INPUT_DIR.exists():
        print(f"错误: 输入目录不存在: {INPUT_DIR}")
        sys.exit(1)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    metrics_list = collect_metrics()
    if RUN_MISSINGDATA_FILTER:
        export_by_thresholds(metrics_list, MISSING_DATA_THRESHOLDS, "missingdata", lambda locus, th: locus.missing_data_ratio <= th, str)
    if RUN_BOOTSTRAP_FILTER:
        export_by_thresholds(metrics_list, BOOTSTRAP_THRESHOLDS, "bootstrap", lambda locus, th: (locus.avg_bootstrap or -1) >= th, str)
    if RUN_SEQ_LENGTH_FILTER:
        export_by_thresholds(metrics_list, SEQ_LENGTH_THRESHOLDS, "seq_length", lambda locus, th: locus.seq_length >= th, str)
    if RUN_PI_SITE_FILTER:
        export_by_thresholds(metrics_list, PI_SITE_THRESHOLDS, "pi_sites", lambda locus, th: locus.pi_sites >= th, str)
    prepare_astral_inputs()
    print(f"筛选完成，输出目录: {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
