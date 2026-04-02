#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INPUT_TREE_FILE = "input.newick"
OUTPUT_TREE_FILE = "support_as_probability.newick"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def convert_support_to_probability(input_file: Path, output_file: Path) -> None:
    newick_string = input_file.read_text(encoding="utf-8")

    def replace_support(match):
        return f"){float(match.group(1)) / 100.0}"

    modified_newick_string = re.sub(r"\)([\d\.]+)", replace_support, newick_string)
    output_file.write_text(modified_newick_string, encoding="utf-8")


def main() -> None:
    input_file = INPUT_DIR / INPUT_TREE_FILE
    if not input_file.exists():
        raise FileNotFoundError(f"未找到输入树文件: {input_file}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    convert_support_to_probability(input_file, OUTPUT_DIR / OUTPUT_TREE_FILE)
    print("树支持率转换完成。")


if __name__ == "__main__":
    main()
