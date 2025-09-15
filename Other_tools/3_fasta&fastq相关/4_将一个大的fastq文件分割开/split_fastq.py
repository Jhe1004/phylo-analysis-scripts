#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
from multiprocessing import Pool, cpu_count

def get_args():
    """
    获取并解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description="一个使用多进程并行分割双端FASTQ文件的脚本 (依赖seqkit)。",
        formatter_class=argparse.RawTextHelpFormatter  # 允许在帮助信息中使用换行符
    )

    # 定义输入文件参数
    parser.add_argument(
        "-p", "--pattern",
        type=str,
        required=True,
        help="用于识别正向测序文件(R1)的模式字符串。\n"
             "例如: 使用 '1.clean.fq.gz' 来匹配像 'SampleA_1.clean.fq.gz' 这样的文件。\n"
             "脚本会自动通过将 '1' 替换为 '2' 来寻找对应的反向测序文件(R2)。"
    )
    
    # 定义进程数参数
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=cpu_count(),  # 默认使用所有可用的 CPU核心数
        help=f"设置并行运行的进程数。\n"
             f"默认值: 当前机器的CPU核心数 ({cpu_count()})。"
    )

    # 定义seqkit内部线程数参数
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=4,
        help="为每个 'seqkit split2' 任务指定的线程数 (即-j参数)。\n"
             "默认值: 4。"
    )
    
    # 定义输出目录前缀
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="split_output",
        help="指定每个分割文件输出目录的前缀。\n"
             "例如: 输入文件为 'SampleA_1.fq.gz'，前缀为 'split_out'，\n"
             "则输出目录将被创建为 'split_out_SampleA'。\n"
             "默认值: 'split_output'。"
    )

    return parser.parse_args()

def get_file_list(pattern):
    """
    获取当前文件夹中符合目标模式的文件。
    
    :param pattern: 用户提供的文件名匹配模式 (e.g., "1.clean.fq.gz")
    :return: file_list: 所有匹配的 R1 文件名列表
    """
    file_list = []
    print(f"正在当前目录中搜索以 '{pattern}' 结尾的文件...")
    for each_file in os.listdir(os.getcwd()):
        if each_file.endswith(pattern):
            # 检查对应的R2文件是否存在
            r2_file = each_file.replace("1.clean.fq.gz", "2.clean.fq.gz")
            if os.path.exists(r2_file):
                file_list.append(each_file)
            else:
                print(f"警告: 找到 R1 文件 '{each_file}' 但未找到对应的 R2 文件 '{r2_file}'。已跳过。")
    
    if not file_list:
        print("错误: 未找到匹配的文件。请检查你的模式字符串和文件位置。")
        exit(1)
        
    print(f"找到 {len(file_list)} 对文件待处理。")
    return file_list

def distribute_tasks(file_list, num_threads):
    """
    将所有文件平均分配给各个进程。

    :param file_list: 所有文件的名称列表
    :param num_threads: 进程数
    :return: task_list: 一个列表，其中每个元素是分配给一个进程的文件列表
    """
    task_list = [[] for _ in range(num_threads)]
    print(f"正在将任务分配到 {num_threads} 个进程中...")
    
    index = 0
    for each_file in file_list:
        task_list[index].append(each_file)
        index = (index + 1) % num_threads
        
    # 过滤掉没有分配到任务的空列表
    task_list = [tasks for tasks in task_list if tasks]
    
    return task_list

def run_split(file_list_chunk, seqkit_jobs, output_prefix):
    """
    对分配到的文件列表执行 seqkit split2 命令。

    :param file_list_chunk: 分配给单个进程的文件名列表
    :param seqkit_jobs: seqkit split2 使用的线程数 (-j)
    :param output_prefix: 输出目录的前缀
    """
    process_id = os.getpid()
    print(f"[进程 {process_id}] 开始处理 {len(file_list_chunk)} 个文件: {', '.join(file_list_chunk)}")
    
    for r1_file in file_list_chunk:
        # 基于R1文件名构建R2文件名和输出目录名
        r2_file = r1_file.replace("1.clean.fq.gz", "2.clean.fq.gz")
        base_name = r1_file.replace("_1.clean.fq.gz", "")
        output_dir = f"{output_prefix}_{base_name}"

        # 构建并执行命令
        command = (
            f"seqkit split2 "
            f"-1 {r1_file} "
            f"-2 {r2_file} "
            f"-p 2 -j {seqkit_jobs} "
            f"-O {output_dir}"
        )
        
        print(f"[进程 {process_id}] 正在执行: {command}")
        try:
            os.system(command)
            print(f"[进程 {process_id}] 成功处理 {r1_file}。")
        except Exception as e:
            print(f"[进程 {process_id}] 处理 {r1_file} 时发生错误: {e}")


def main():
    """
    主函数
    """
    # 1. 获取参数
    args = get_args()
    
    # 2. 获取文件列表
    file_list = get_file_list(args.pattern)
    
    # 3. 如果文件数小于进程数，调整进程数以避免资源浪费
    num_threads = min(args.threads, len(file_list))
    if num_threads == 0:
        print("没有文件需要处理，程序退出。")
        return
        
    # 4. 分配任务
    task_list = distribute_tasks(file_list, num_threads)
    
    # 5. 创建并运行进程池
    # 使用 with 语句可以确保进程池在使用后被正确关闭
    with Pool(processes=num_threads) as p:
        # 使用 starmap 可以方便地传递多个参数
        tasks_to_run = [(chunk, args.jobs, args.output_prefix) for chunk in task_list]
        p.starmap(run_split, tasks_to_run)

    print("\n----- 所有任务已完成。 -----")


if __name__ == "__main__":
    main()