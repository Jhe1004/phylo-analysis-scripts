import rasterio
from rasterio.warp import reproject, Resampling
from pathlib import Path
import os
import sys

# --- 1. 配置您的路径 ---
try:
    DATA_DIR = Path(__file__).parent
except NameError:
    DATA_DIR = Path.cwd()
    print("⚠️ 警告：在交互模式下运行，使用当前工作目录。")
    print(f"   请确保 {DATA_DIR} 包含了您的数据文件。")

OUTPUT_DIR_NAME = "output_masked"
REF_FILE_NAME = "ref.tif"
SOURCE_EXTENSIONS = [".asc"]

# --- 2. 脚本正文 ---

def process_rasters(data_dir, output_dir_name, ref_file_name, source_extensions):
    """
    使用参考栅格文件对文件夹中的所有源栅格文件进行掩膜提取和重采样。
    """
    data_path = Path(data_dir)
    ref_file_path = data_path / ref_file_name
    output_path = data_path / output_dir_name

    print(f"ℹ️ 脚本启动。")
    print(f"  > 数据目录: {data_path.resolve()}")

    try:
        output_path.mkdir(parents=True, exist_ok=True)
        print(f"✅ 输出文件夹已创建/确认: {output_path.resolve()}")
    except Exception as e:
        print(f"❌ 无法创建输出文件夹: {e}")
        return

    if not ref_file_path.exists():
        print(f"❌ 错误：找不到参考文件 {ref_file_path}")
        print(f"  > 请确保 '{ref_file_name}' 与 .py 脚本在同一个文件夹中。")
        return

    print(f"\n[1/3] 正在加载参考图层: {ref_file_name}...")
    try:
        with rasterio.open(ref_file_path) as ref_ds:
            dst_crs = ref_ds.crs
            dst_transform = ref_ds.transform
            dst_width = ref_ds.width
            dst_height = ref_ds.height
            dst_nodata = ref_ds.nodata
            
            # ✨ --- (修复 1/1) --- 
            # 从 .dtype (错误) 更改为 .dtypes[0] (正确)
            dst_dtype = ref_ds.dtypes[0] 
            
            print(f"  > 参考图层属性加载成功。")
            print(f"  > 目标 CRS: {dst_crs}")
            print(f"  > 目标尺寸: ({dst_width} W, {dst_height} H)")
            print(f"  > 目标 Dtype: {dst_dtype}") # 打印出来确认

    except Exception as e:
        print(f"❌ 打开参考文件 {ref_file_name} 时出错: {e}")
        return

    print(f"\n[2/3] 正在查找源文件 (扩展名: {', '.join(source_extensions)})...")
    all_source_files = []
    for ext in source_extensions:
        all_source_files.extend(data_path.glob(f"*{ext}"))

    ref_path_resolved = ref_file_path.resolve()
    files_to_process = [
        f for f in all_source_files if f.resolve() != ref_path_resolved
    ]

    if not files_to_process:
        print(f"⚠️ 警告：在 {data_path} 中未找到任何匹配的源文件。")
        return

    print(f"  > 找到 {len(files_to_process)} 个文件准备处理。")

    count_success = 0
    count_fail = 0
    for src_file in files_to_process:
        try:
            print(f"  > 正在处理: {src_file.name}...")
            with rasterio.open(src_file) as src_ds:
                
                out_meta = src_ds.meta.copy()
                
                out_meta.update({
                    "driver": "GTiff",
                    "crs": dst_crs,
                    "transform": dst_transform,
                    "width": dst_width,
                    "height": dst_height,
                    "nodata": dst_nodata,
                    "dtype": dst_dtype # 这里使用上面获取的 dst_dtype
                })

                output_file_path = output_path / (src_file.stem + '_masked.tif')

                with rasterio.open(output_file_path, 'w', **out_meta) as dst_ds:
                    reproject(
                        source=rasterio.band(src_ds, 1),
                        destination=rasterio.band(dst_ds, 1),
                        src_transform=src_ds.transform,
                        src_crs=src_ds.crs,
                        dst_transform=dst_transform,
                        dst_crs=dst_crs,
                        resampling=Resampling.bilinear
                    )
            print(f"    ... 成功保存到: {output_file_path.name}")
            count_success += 1

        except Exception as e:
            print(f"    ... ❌ 处理 {src_file.name} 时发生错误: {e}")
            count_fail += 1

    print(f"\n[3/3] ✨ 处理完成！")
    print(f"  > {count_success} 个文件成功。")
    print(f"  > {count_fail} 个文件失败。")
    print(f"  > 所有文件已保存至: {output_path.resolve()}")

# --- 3. 运行脚本 ---
if __name__ == "__main__":
    if 'rasterio' not in sys.modules:
        print("错误：未找到 'rasterio' 库。")
        print("请使用 'conda install -c conda-forge rasterio' 或 'pip install rasterio' 安装。")
    else:
        process_rasters(DATA_DIR, OUTPUT_DIR_NAME, REF_FILE_NAME, SOURCE_EXTENSIONS)