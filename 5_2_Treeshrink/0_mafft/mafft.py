import os
import multiprocessing

# 定义全局变量
file_tag = ".fasta"  # 目标文件的扩展名
num_cpu = 1  # 并行进程数（一次处理的文件数）

def get_file_list():
    """
    获取当前文件夹中符合目标扩展名的文件列表
    """
    file_name = []
    for each in os.listdir(os.getcwd()):
        if file_tag in each:
            file_name.append(each)
    return file_name

def run_mafft(each_file):
    """
    对单个文件运行MAFFT
    """
    thread = 15  # 每个进程的线程数，可以根据需要调整
    output_file = each_file.replace(file_tag, "") + "_maffted.fasta"
    command = (f"mafft --adjustdirection --auto --thread {thread} {each_file} | sed 's/_R_//g' > {output_file}")
    
    try:
        print(f"正在处理文件: {each_file}")
        os.system(command)
        print(f"处理完成: {each_file}")
    except Exception as e:
        print(f"处理文件 {each_file} 时出错: {e}")

def main():
    # 获取目标文件列表
    file_list = get_file_list()
    
    # 创建进程池并并行处理文件
    with multiprocessing.Pool(processes=num_cpu) as pool:
        pool.map(run_mafft, file_list)

if __name__ == "__main__":
    main()
