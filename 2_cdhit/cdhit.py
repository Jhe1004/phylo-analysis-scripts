import os
import subprocess
from multiprocessing import Pool, cpu_count
import logging

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

def get_file_list():
    """
    获取当前目录下所有以 .fas 结尾的文件名。
    """
    now_dir = os.getcwd()
    file_temp = os.listdir(path=now_dir)
    file_name = [file for file in file_temp if file.endswith(".fas")]
    return file_name

def process_fasta(fasta_name):
    """
    对单个 .fas 文件执行 cd-hit-est 命令。
    """
    try:
        input_file = fasta_name
        # 构建输出文件名，避免在文件名中出现多余的空格
        output_file = f"{fasta_name}ta"

        # 构建命令
        cmd = ["./cd-hit-est", "-i", input_file, "-o", output_file]

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

    # 获取所有 .fas 文件
    file_names = get_file_list()
    if not file_names:
        logging.warning("No .fas files found in the current directory.")
        return

    logging.info(f"Found {len(file_names)} .fas files to process.")

    # 设置进程池的大小（默认为 CPU 核心数）
    pool_size = cpu_count()
    logging.info(f"Starting processing with {pool_size} parallel processes...")

    # 创建进程池并并行处理文件
    with Pool(pool_size) as pool:
        pool.map(process_fasta, file_names)

    logging.info("All files have been processed.")

if __name__ == "__main__":
    main()