#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hashlib
from pathlib import Path


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INPUT_FILE_LIST = "md5.txt"
OUTPUT_TSV_FILE = "md5_results.tsv"
BLOCK_SIZE = 65536


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def get_file_md5(filename: Path, block_size: int = BLOCK_SIZE) -> str | None:
    md5_hash = hashlib.md5()
    try:
        with filename.open("rb") as f:
            for chunk in iter(lambda: f.read(block_size), b""):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()
    except IOError:
        return None


def main() -> None:
    input_list_file = INPUT_DIR / INPUT_FILE_LIST
    if not input_list_file.exists():
        raise FileNotFoundError(f"未找到文件列表: {input_list_file}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    filenames = [line.strip() for line in input_list_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    output_tsv_file = OUTPUT_DIR / OUTPUT_TSV_FILE
    with output_tsv_file.open("w", encoding="utf-8") as out:
        for filename in filenames:
            file_path = INPUT_DIR / filename
            if not file_path.is_file():
                md5_sum = "FILE_NOT_FOUND"
            else:
                md5_sum = get_file_md5(file_path) or "READ_ERROR"
            out.write(f"{filename}\t{md5_sum}\n")
    print("MD5 计算完成。")


if __name__ == "__main__":
    main()
