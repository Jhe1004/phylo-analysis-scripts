import os
from multiprocessing import Pool, cpu_count

def get_file_list():
    """
    获取当前目录下所有以 .fas 结尾的文件名。
    """
    now_dir = os.getcwd()           # 获取当前工作目录
    file_temp = os.listdir(path=now_dir)    # 列出当前目录中的所有文件
    file_name = [file for file in file_temp if file.endswith(".fasta")]  # 筛选出 .fas 文件
    return file_name

def process_fasta(fasta_name):
    """
    对单个 .fas 文件执行 TransDecoder.LongOrfs 和 TransDecoder.Predict 命令。
    """
    try:
        # 执行 TransDecoder.LongOrfs 命令
        cmd_longorfs = f"./TransDecoder.LongOrfs -t {fasta_name}"
        os.system(cmd_longorfs)
        print(f"Completed LongOrfs for {fasta_name}")

        # 执行 TransDecoder.Predict 命令
        cmd_predict = f"./TransDecoder.Predict --no_refine_starts -t {fasta_name}"
        os.system(cmd_predict)
        print(f"Completed Predict for {fasta_name}")

    except Exception as e:
        print(f"Error processing {fasta_name}: {e}")

def main():
    # 获取所有 .fas 文件
    file_names = get_file_list()
    if not file_names:
        print("No .fas files found in the current directory.")
        return

    # 设置进程池的大小（默认为 CPU 核心数）
    pool_size = cpu_count()
    print(f"Starting processing with {pool_size} parallel processes...")

    # 创建进程池并并行处理文件
    with Pool(pool_size) as pool:
        pool.map(process_fasta, file_names)

    print("All files have been processed.")

if __name__ == "__main__":
    main()