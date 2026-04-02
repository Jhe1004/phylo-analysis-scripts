#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


READ1_SUFFIXES = (
    "_1.sub10.fq.gz",
    "_1.clean.fq.gz",
    "_1.fastq.gz",
    "_1.fq.gz",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="批量调用 GetOrganelle 组装 nrDNA（双端测序）"
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        default=".",
        help="输入目录，默认当前目录",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="getorganelle_nrdna_results",
        help="输出目录，默认: getorganelle_nrdna_results",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="每个样本使用的线程数，默认: 16",
    )
    parser.add_argument(
        "--env-name",
        default="getorganelle",
        help="conda 环境名，默认: getorganelle",
    )
    parser.add_argument(
        "--getorg-entry",
        default=None,
        help="GetOrganelle 程序名，可选 get_organelle_from_reads.py 或 get_organelle_from_reads",
    )
    parser.add_argument(
        "-F",
        "--organelle-type",
        default="embplant_nr",
        help="GetOrganelle 的 -F 参数，默认: embplant_nr",
    )
    parser.add_argument(
        "-k",
        "--kmer-list",
        default="35,85,115",
        help="GetOrganelle 的 -k 参数，默认: 35,85,115",
    )
    parser.add_argument(
        "-R",
        "--max-rounds",
        type=int,
        default=10,
        help="GetOrganelle 的 -R 参数，默认: 10",
    )
    parser.add_argument(
        "--extra-args",
        default="",
        help="额外传给 GetOrganelle 的参数，例如 \"-s seed.fa --max-reads 1.5e8\"",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="即使检测到已有结果，也强制重新运行",
    )
    return parser.parse_args()


def detect_entry(env_name: str, user_entry: str | None) -> str:
    candidates = [user_entry] if user_entry else [
        "get_organelle_from_reads.py",
        "get_organelle_from_reads",
    ]
    for entry in candidates:
        if not entry:
            continue
        cmd = ["conda", "run", "-n", env_name, entry, "-h"]
        result = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if result.returncode == 0:
            return entry
    raise RuntimeError(f"在 conda 环境 {env_name} 中没有找到 GetOrganelle 主程序。")


def find_pairs(input_dir: Path) -> list[tuple[str, Path, Path]]:
    pairs: list[tuple[str, Path, Path]] = []
    seen = set()
    for path in sorted(input_dir.iterdir(), key=lambda p: p.name):
        if not path.is_file():
            continue
        name = path.name
        matched_suffix = next((s for s in READ1_SUFFIXES if name.endswith(s)), None)
        if not matched_suffix:
            continue
        sample = name[: -len(matched_suffix)]
        mate_name = f"{sample}{matched_suffix.replace('_1', '_2', 1)}"
        mate = input_dir / mate_name
        if mate.is_file() and sample not in seen:
            pairs.append((sample, path, mate))
            seen.add(sample)
    return pairs


def result_exists(sample_out: Path) -> bool:
    return (
        (sample_out / "extended_path_sequence.fasta").is_file()
        or (sample_out / "final_assembly_graph.fastg").is_file()
    )


def main() -> int:
    args = parse_args()

    if shutil.which("conda") is None:
        print("错误：未找到 conda，请先确保 conda 可用。", file=sys.stderr)
        return 1

    input_dir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()
    logdir = outdir / "logs"
    outdir.mkdir(parents=True, exist_ok=True)
    logdir.mkdir(parents=True, exist_ok=True)

    try:
        entry = detect_entry(args.env_name, args.getorg_entry)
    except RuntimeError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    pairs = find_pairs(input_dir)
    if not pairs:
        print("当前目录未找到任何 *_1.fq.gz / *_1.fastq.gz 类型的双端输入文件。", file=sys.stderr)
        return 1

    extra_args = shlex.split(args.extra_args)

    print(f"输入目录: {input_dir}")
    print(f"输出目录: {outdir}")
    print(f"日志目录: {logdir}")
    print(f"线程数: {args.threads}")
    print(f"GetOrganelle 命令: conda run -n {args.env_name} {entry}")
    print(f"组装类型: {args.organelle_type}")
    print(f"k-mer: {args.kmer_list}")
    print(f"最大迭代轮数: {args.max_rounds}")
    print()

    done_n = 0
    skip_n = 0
    fail_n = 0

    for sample, fq1, fq2 in pairs:
        sample_out = outdir / sample
        sample_log = logdir / f"{sample}.log"

        if not args.force and result_exists(sample_out):
            print(f"已存在结果，跳过: {sample}")
            skip_n += 1
            continue

        sample_out.mkdir(parents=True, exist_ok=True)

        cmd = [
            "conda",
            "run",
            "-n",
            args.env_name,
            entry,
            "-1",
            str(fq1),
            "-2",
            str(fq2),
            "-o",
            str(sample_out),
            "-F",
            args.organelle_type,
            "-R",
            str(args.max_rounds),
            "-k",
            args.kmer_list,
            "-t",
            str(args.threads),
            *extra_args,
        ]

        print("=" * 60)
        print(f"开始样本: {sample}")
        print(f"R1: {fq1}")
        print(f"R2: {fq2}")
        print(f"输出: {sample_out}")
        print(f"日志: {sample_log}")

        with sample_log.open("w", encoding="utf-8") as log_handle:
            log_handle.write("命令:\n")
            log_handle.write(" ".join(shlex.quote(x) for x in cmd))
            log_handle.write("\n\n")
            result = subprocess.run(
                cmd,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                check=False,
            )

        if result.returncode == 0:
            print(f"完成: {sample}")
            done_n += 1
        else:
            print(f"失败: {sample}，请检查日志 {sample_log}")
            fail_n += 1

    print()
    print("全部任务结束")
    print(f"成功: {done_n}")
    print(f"跳过: {skip_n}")
    print(f"失败: {fail_n}")

    return 1 if fail_n else 0


if __name__ == "__main__":
    sys.exit(main())
