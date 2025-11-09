import rasterio
from pathlib import Path
import os
import sys

# --- 1. 请在这里配置您的转换任务 ---

# 自动获取当前.py脚本所在的文件夹路径
try:
    DATA_DIR = Path(__file__).parent
except NameError:
    DATA_DIR = Path.cwd()
    print("⚠️ 警告：在交互模式下运行，使用当前工作目录。")

# 您希望保存转换后图层的文件夹名称
OUTPUT_DIR_NAME = "output_converted"

# (A) 您想要转换的文件的扩展名 (可以有多个)
# 示例: [".tif", ".img"]
INPUT_EXTENSIONS = [".tif"]

# (B) 您想要转换成的目标格式 (一次只能选一个)
# 请从下面的 'KEY' 中选择一个:
#   'TIF'  -> GeoTIFF (.tif)
#   'ASC'  -> Esri ASCII Grid (.asc)
#   'BIL'  -> Esri BIL (.bil + .hdr)
#   'IMG'  -> ERDAS Imagine (.img)
#
OUTPUT_FORMAT_KEY = "ASC"

# ----------------------------------------------------
# --- 2. 驱动映射 (一般无需修改) ---
# ----------------------------------------------------

# 将用户友好的 KEY 映射到 rasterio/GDAL 内部的驱动名称和新扩展名
FORMAT_MAP = {
    "TIF":  {"driver": "GTiff",   "ext": ".tif"},
    "ASC":  {"driver": "AAIGrid", "ext": ".asc"},
    "BIL":  {"driver": "EHdr",    "ext": ".bil"},
    "IMG":  {"driver": "HFA",     "ext": ".img"},
    # 您可以根据 GDAL 驱动列表按需添加更多
}

# ----------------------------------------------------
# --- 3. 脚本正文 ---
# ----------------------------------------------------

def convert_rasters(data_dir, output_dir_name, input_exts, output_key):
    """
    在同一文件夹内批量转换栅格文件格式。
    """
    data_path = Path(data_dir)
    output_path = data_path / output_dir_name

    print(f"ℹ️ 脚本启动。")
    print(f"  > 数据目录: {data_path.resolve()}")

    # 1. 验证输出格式
    if output_key not in FORMAT_MAP:
        print(f"❌ 错误：输出格式 '{output_key}' 不被支持。")
        print(f"  > 请从 {list(FORMAT_MAP.keys())} 中选择。")
        return
    
    driver_config = FORMAT_MAP[output_key]
    out_driver = driver_config["driver"]
    out_ext = driver_config["ext"]
    
    print(f"  > 转换任务: {', '.join(input_exts)} -> {output_key} ({out_ext})")

    # 2. 创建输出文件夹
    try:
        output_path.mkdir(parents=True, exist_ok=True)
        print(f"✅ 输出文件夹已创建/确认: {output_path.resolve()}")
    except Exception as e:
        print(f"❌ 无法创建输出文件夹: {e}")
        return

    # 3. 查找所有要处理的源文件
    files_to_process = []
    for ext in input_exts:
        files_to_process.extend(data_path.glob(f"*{ext}"))

    if not files_to_process:
        print(f"⚠️ 警告：在 {data_path} 中未找到任何匹配 '{' '.join(input_exts)}' 的文件。")
        return

    print(f"\n[1/2] 找到 {len(files_to_process)} 个文件准备转换...")

    # 4. 循环处理每个文件
    count_success = 0
    count_fail = 0
    for src_file in files_to_process:
        print(f"  > 正在转换: {src_file.name}...")
        try:
            with rasterio.open(src_file) as src_ds:
                # 读取源文件的元数据和数据
                meta = src_ds.meta.copy()
                data = src_ds.read() # 读取所有波段

                # 检查特殊格式的限制
                if out_driver == 'AAIGrid' and src_ds.count > 1:
                    print(f"    ... ⚠️ 跳过: {src_file.name} 是多波段文件。")
                    print(f"    ... 'ASC' 格式 (AAIGrid) 只支持单波段文件。")
                    count_fail += 1
                    continue

                # 更新元数据以使用新的驱动
                meta.update(driver=out_driver)
                
                # 定义新的输出文件名
                output_file_path = output_path / (src_file.stem + out_ext)

                # 写入新文件
                with rasterio.open(output_file_path, 'w', **meta) as dst_ds:
                    dst_ds.write(data)
            
            print(f"    ... 成功保存到: {output_file_path.name}")
            count_success += 1

        except Exception as e:
            print(f"    ... ❌ 处理 {src_file.name} 时发生错误: {e}")
            count_fail += 1

    print(f"\n[2/2] ✨ 转换完成！")
    print(f"  > {count_success} 个文件成功。")
    print(f"  > {count_fail} 个文件失败或被跳过。")
    print(f"  > 所有文件已保存至: {output_path.resolve()}")

# --- 4. 运行脚本 ---
if __name__ == "__main__":
    if 'rasterio' not in sys.modules:
        print("错误：未找到 'rasterio' 库。")
        print("请使用 'conda install -c conda-forge rasterio' 或 'pip install rasterio' 安装。")
    else:
        convert_rasters(DATA_DIR, OUTPUT_DIR_NAME, INPUT_EXTENSIONS, OUTPUT_FORMAT_KEY)