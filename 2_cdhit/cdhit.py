import os
import subprocess
import shutil
from multiprocessing import Pool, cpu_count
import logging
import sys


def setup_logging():
    """
    设置日志记录配置。
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler("processing.log"),
            logging.StreamHandler()
        ]
    )

def find_cdhit_executable():
    """
    查找 cd-hit-est 可执行文件的路径。
    """
    # 首先检查命令行参数是否指定了路径
    if len(sys.argv) > 1:
        cdhit_path = sys.argv[1]
        if os.path.isfile(cdhit_path) and os.access(cdhit_path, os.X_OK):
            return cdhit_path
        else:
            logging.error(f"指定的 cd-hit-est 路径无效或不可执行: {cdhit_path}")
            return None
    
    # 检查系统路径中是否有 cd-hit-est
    cdhit_path = shutil.which("cd-hit-est")
    if cdhit_path is not None:
        return cdhit_path
    
    # 检查当前目录是否有 cd-hit-est
    if os.path.isfile("./cd-hit-est") and os.access("./cd-hit-est", os.X_OK):
        return "./cd-hit-est"
    
    logging.error("cd-hit-est 未安装或不在系统路径中。")
    logging.error("请先安装 cd-hit 软件包，或确保 cd-hit-est 可执行文件在系统路径中。")
    logging.error("或者作为参数提供 cd-hit-est 的完整路径，例如: python cdhit.py /path/to/cd-hit-est")
    return None

def get_file_list():
    """
    获取当前目录下所有以 .fas 结尾的文件名。
    """
    now_dir = os.getcwd()
    file_temp = os.listdir(path=now_dir)
    file_name = [file for file in file_temp if file.endswith(".fas")]
    return file_name

def process_fasta(args):
    """
    对单个 .fas 文件执行 cd-hit-est 命令。
    """
    fasta_name, cdhit_path = args
    try:
        input_file = fasta_name
        # 构建输出文件名，避免在文件名中出现多余的空格
        output_file = f"{fasta_name}ta"

        # 构建命令
        cmd = [cdhit_path, "-i", input_file, "-o", output_file]

        logging.info(f"Starting cd-hit-est for {input_file} -> {output_file}")

        # 执行命令
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # 检查命令是否成功
        if result.returncode == 0:
            logging.info(f"Successfully processed {input_file}")
        else:
            logging.error(f"Error processing {input_file}: {result.stderr}")

    except Exception as e:
        logging.error(f"Exception processing {fasta_name}: {e}")

def main():
    setup_logging()

    # 查找 cd-hit-est 可执行文件
    cdhit_path = find_cdhit_executable()
    if cdhit_path is None:
        return

    # 获取所有 .fas 文件
    file_names = get_file_list()
    if not file_names:
        logging.warning("No .fas files found in the current directory.")
        return

    logging.info(f"Found {len(file_names)} .fas files to process.")

    # 设置进程池的大小（默认为 CPU 核心数）
    pool_size = cpu_count()
    logging.info(f"Starting processing with {pool_size} parallel processes...")

    # 创建参数列表，包含文件名和 cd-hit-est 路径
    args_list = [(file_name, cdhit_path) for file_name in file_names]

    # 创建进程池并并行处理文件
    with Pool(pool_size) as pool:
        pool.map(process_fasta, args_list)

    logging.info("All files have been processed.")

if __name__ == "__main__":
    main()