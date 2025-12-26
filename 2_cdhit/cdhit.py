import os
import subprocess
import shutil
from multiprocessing import Pool, cpu_count
import logging
import sys

# --- 配置参数 ---
SIF_FILENAME = "cdhit.sif"  # 镜像文件名
CDHIT_CMD = "cd-hit-est"    # 容器内要调用的命令 (如果是蛋白序列请改为 cd-hit)
IDENTITY_THRESHOLD = "0.95" # 相似度阈值 (-c 参数)
# ----------------

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

def find_sif_image():
    """
    查找 cdhit.sif 镜像文件的路径。
    """
    if len(sys.argv) > 1:
        sif_path = sys.argv[1]
        if os.path.isfile(sif_path) and sif_path.endswith(".sif"):
            return os.path.abspath(sif_path)

    current_dir_sif = os.path.join(os.getcwd(), SIF_FILENAME)
    if os.path.isfile(current_dir_sif):
        return current_dir_sif
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir_sif = os.path.join(script_dir, SIF_FILENAME)
    if os.path.isfile(script_dir_sif):
        return script_dir_sif

    logging.error(f"未找到镜像文件: {SIF_FILENAME}")
    logging.error(f"请确保 {SIF_FILENAME} 在当前目录，或者通过参数指定完整路径。")
    return None

def check_apptainer_installed():
    """
    检查系统是否安装了 apptainer。
    """
    if shutil.which("apptainer") is None:
        logging.error("未找到 'apptainer' 命令。请先安装 Apptainer/Singularity。")
        return False
    return True

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
    使用 Apptainer 镜像对单个 .fas 文件执行 cd-hit 命令。
    """
    fasta_name, sif_path = args
    try:
        input_file = fasta_name
        # 构建输出文件名
        output_file = f"{fasta_name}ta" 

        # --- 构建 Apptainer 命令 ---
        cmd = [
            "apptainer", "exec", 
            sif_path, 
            CDHIT_CMD, 
            "-i", input_file, 
            "-o", output_file, 
            "-c", IDENTITY_THRESHOLD
        ]

        logging.info(f"Starting: {' '.join(cmd)}")

        # 执行命令
        # 修正点：将 text=True 改为 universal_newlines=True 以兼容 Python 3.6
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

        if result.returncode == 0:
            logging.info(f"Successfully processed {input_file}")
        else:
            logging.error(f"Error processing {input_file}:")
            logging.error(f"STDOUT: {result.stdout}")
            logging.error(f"STDERR: {result.stderr}")

    except Exception as e:
        logging.error(f"Exception processing {fasta_name}: {e}")

def main():
    setup_logging()

    if not check_apptainer_installed():
        return

    sif_path = find_sif_image()
    if sif_path is None:
        return
    logging.info(f"Using Apptainer image: {sif_path}")

    file_names = get_file_list()
    if not file_names:
        logging.warning("No .fas files found in the current directory.")
        return

    logging.info(f"Found {len(file_names)} .fas files to process.")

    pool_size = cpu_count()
    logging.info(f"Starting processing with {pool_size} parallel processes...")

    args_list = [(file_name, sif_path) for file_name in file_names]

    with Pool(pool_size) as pool:
        pool.map(process_fasta, args_list)

    logging.info("All files have been processed.")

if __name__ == "__main__":
    main()