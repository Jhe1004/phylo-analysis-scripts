#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import re
import shutil


SCRIPT_DIR = Path(__file__).resolve().parent

# ==============================
# 配置区：请直接修改这里
# ==============================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"

INPUT_EXTENSION = ".gb"
COPY_UNMODIFIED_FILES = True
# ==============================


def main() -> None:
    input_dir = SCRIPT_DIR / INPUT_DIRECTORY
    output_dir = SCRIPT_DIR / OUTPUT_DIRECTORY
    output_dir.mkdir(parents=True, exist_ok=True)

    gb_files = sorted(input_dir.glob(f"*{INPUT_EXTENSION}"))
    if not gb_files:
        raise FileNotFoundError(f"在 {input_dir} 中未找到 {INPUT_EXTENSION} 文件")

    for gb_file in gb_files:
        new_locus_name = re.sub(r"_+", "_", re.sub(r"[^a-zA-Z0-9]", "_", gb_file.stem)).strip("_")
        lines = gb_file.read_text(encoding="utf-8").splitlines(keepends=True)
        if not lines:
            continue
        if lines[0].startswith("LOCUS"):
            parts = lines[0].split()
            if len(parts) > 1:
                parts[1] = new_locus_name
                lines[0] = " ".join(parts) + "\n"
        output_file = output_dir / gb_file.name
        output_file.write_text("".join(lines), encoding="utf-8")

    if COPY_UNMODIFIED_FILES:
        for other_file in input_dir.iterdir():
            if other_file.is_file() and other_file.suffix != INPUT_EXTENSION:
                shutil.copy2(other_file, output_dir / other_file.name)

    print(f"完成。输出目录: {output_dir}")


if __name__ == "__main__":
    main()
