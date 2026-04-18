#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

import rasterio
from rasterio.warp import Resampling, reproject


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/4_arcgis相关/0_掩膜提取/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/4_arcgis相关/0_掩膜提取/output"
REF_FILE_NAME = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/4_arcgis相关/0_掩膜提取/input/ref.tif"
SOURCE_EXTENSIONS = [".asc"]


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    ref_file_path = INPUT_DIR / REF_FILE_NAME
    if not ref_file_path.exists():
        raise FileNotFoundError(f"未找到参考图层: {ref_file_path}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with rasterio.open(ref_file_path) as ref_ds:
        dst_crs = ref_ds.crs
        dst_transform = ref_ds.transform
        dst_width = ref_ds.width
        dst_height = ref_ds.height
        dst_nodata = ref_ds.nodata
        dst_dtype = ref_ds.dtypes[0]
    files_to_process = []
    for ext in SOURCE_EXTENSIONS:
        files_to_process.extend(sorted(INPUT_DIR.glob(f"*{ext}")))
    files_to_process = [f for f in files_to_process if f.resolve() != ref_file_path.resolve()]
    for src_file in files_to_process:
        with rasterio.open(src_file) as src_ds:
            out_meta = src_ds.meta.copy()
            out_meta.update(driver="GTiff", crs=dst_crs, transform=dst_transform, width=dst_width, height=dst_height, nodata=dst_nodata, dtype=dst_dtype)
            output_file_path = OUTPUT_DIR / f"{src_file.stem}_masked.tif"
            with rasterio.open(output_file_path, "w", **out_meta) as dst_ds:
                reproject(
                    source=rasterio.band(src_ds, 1),
                    destination=rasterio.band(dst_ds, 1),
                    src_transform=src_ds.transform,
                    src_crs=src_ds.crs,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.bilinear,
                )
    print("掩膜提取完成。")


if __name__ == "__main__":
    main()
