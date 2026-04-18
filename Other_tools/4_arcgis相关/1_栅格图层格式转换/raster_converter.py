#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

import rasterio


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/4_arcgis相关/1_栅格图层格式转换/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/4_arcgis相关/1_栅格图层格式转换/output"
INPUT_EXTENSIONS = [".tif"]
OUTPUT_FORMAT_KEY = "ASC"

FORMAT_MAP = {
    "TIF": {"driver": "GTiff", "ext": ".tif"},
    "ASC": {"driver": "AAIGrid", "ext": ".asc"},
    "BIL": {"driver": "EHdr", "ext": ".bil"},
    "IMG": {"driver": "HFA", "ext": ".img"},
}


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    if OUTPUT_FORMAT_KEY not in FORMAT_MAP:
        raise ValueError(f"不支持的输出格式: {OUTPUT_FORMAT_KEY}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_driver = FORMAT_MAP[OUTPUT_FORMAT_KEY]["driver"]
    out_ext = FORMAT_MAP[OUTPUT_FORMAT_KEY]["ext"]
    files_to_process = []
    for ext in INPUT_EXTENSIONS:
        files_to_process.extend(sorted(INPUT_DIR.glob(f"*{ext}")))
    if not files_to_process:
        raise FileNotFoundError(f"未在 {INPUT_DIR} 中找到待转换栅格文件。")
    for src_file in files_to_process:
        with rasterio.open(src_file) as src_ds:
            meta = src_ds.meta.copy()
            data = src_ds.read()
            if out_driver == "AAIGrid" and src_ds.count > 1:
                continue
            meta.update(driver=out_driver)
            with rasterio.open(OUTPUT_DIR / f"{src_file.stem}{out_ext}", "w", **meta) as dst_ds:
                dst_ds.write(data)
    print("栅格格式转换完成。")


if __name__ == "__main__":
    main()
