#!/usr/bin/env python3
import os
import subprocess
import argparse
import sys
# 导入 multiprocessing 模块中的 Pool 和 cpu_count
from multiprocessing import Pool, cpu_count

def find_read_pairs(directory, suffix1, suffix2):
    """
    在指定目录中查找成对的测序文件。
    （此函数未作修改）
    """
    pairs = []
    for filename in os.listdir(directory):
        if filename.endswith(suffix1):
            prefix = filename[:-len(suffix1)]
            read1_path = os.path.join(directory, filename)
            read2_filename = prefix + suffix2
            read2_path = os.path.join(directory, read2_filename)

            if os.path.exists(read2_path):
                pairs.append((prefix, read1_path, read2_path))
            else:
                print(f"⚠️  警告: 找到了文件 '{filename}' 但未找到其配对文件 '{read2_filename}'。将跳过此样本。")
    
    if not pairs:
        print(f"❌ 错误: 在当前目录下未找到任何以 '{suffix1}' 结尾的文件。请检查后缀名是否正确。")
        sys.exit(1)
        
    return pairs

def run_getorganelle(sample_prefix, read1, read2, assembly_type, threads):
    """
    为单个样本构建并执行 GetOrganelle 命令。
    （此函数未作修改，它将由每个子进程独立调用）
    """
    if assembly_type == 'plastome':
        output_dir = f"{sample_prefix}_plastome_out"
        organelle_type = "embplant_pt"
        rounds = "30"
        kmers = "21,45,65,85,105,125"
        print(f"\n assembling plastome for {sample_prefix}...")
    elif assembly_type == 'rdna':
        output_dir = f"{sample_prefix}_nrDNA_out"
        organelle_type = "embplant_nr"
        rounds = "10"
        kmers = "35,85,115"
        print(f"\n assembling nrDNA for {sample_prefix}...")
    else:
        raise ValueError("无效的组装类型。")

    command = [
        "get_organelle_from_reads.py",
        "-1", read1,
        "-2", read2,
        "-o", output_dir,
        "-R", rounds,
        "-k", kmers,
        "-F", organelle_type,
        "-t", str(threads),
    ]

    # 打印信息加上进程ID，方便区分不同进程的输出
    pid = os.getpid()
    print("-" * 60)
    print(f" [Process PID: {pid}] Executing command for sample: {sample_prefix}")
    print(f" {' '.join(command)}")
    print("-" * 60)

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, encoding='utf-8')
        print(f"--- STDOUT from {sample_prefix} (PID: {pid}) ---\n{result.stdout}")
        print(f"✅ 成功 (PID: {pid}): {sample_prefix} 的组装任务已完成。结果保存在 '{output_dir}' 目录中。")
    except FileNotFoundError:
        print(f"❌ 错误 (PID: {pid}): 'get_organelle_from_reads.py' 命令未找到。")
        print("请确保 GetOrganelle 已正确安装并且其路径已添加到系统的 PATH 环境变量中。")
    except subprocess.CalledProcessError as e:
        print(f"❌ 错误 (PID: {pid}): {sample_prefix} 的组装任务失败。")
        print(f"Return code: {e.returncode}")
        print(f"\n--- GetOrganelle STDOUT (PID: {pid}) ---")
        print(e.stdout)
        print(f"\n--- GetOrganelle STDERR (PID: {pid}) ---")
        print(e.stderr)
        print("---------------------------\n")


def main():
    """主函数，解析参数并以多进程方式启动流程。"""
    # --- MODIFIED: 增加了 -p/--processes 参数 ---
    default_processes = min(4, cpu_count()) # 默认使用4个进程或CPU核心数中较小的值

    parser = argparse.ArgumentParser(
        description="以多进程方式批量运行 GetOrganelle 进行叶绿体或 nrDNA 组装。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "-t", "--type",
        required=True,
        choices=['plastome', 'rdna'],
        help="必须指定。选择要组装的基因组类型:\n"
             "  'plastome' - 组装叶绿体基因组 (使用 embplant_pt)。\n"
             "  'rdna'     - 组装核糖体DNA (使用 embplant_nr)。"
    )
    parser.add_argument(
        "-s1", "--suffix1",
        required=True,
        help="必须指定。第一个 (Forward/R1) 测序文件的后缀。\n"
             "例如: '_1.clean.fq.gz' 或 '_R1.fastq'"
    )
    parser.add_argument(
        "-s2", "--suffix2",
        required=True,
        help="必须指定。第二个 (Reverse/R2) 测序文件的后缀。\n"
             "例如: '_2.clean.fq.gz' 或 '_R2.fastq'"
    )
    
    parser.add_argument(
        "-T", "--threads",
        type=int,
        default=4,
        help="指定每个GetOrganelle任务运行所使用的线程数 (传递给其-t参数)。\n"
             "默认值: 4。"
    )
    # --- NEW: 添加 -p/--processes 参数用于控制并行进程数 ---
    parser.add_argument(
        "-p", "--processes",
        type=int,
        default=default_processes,
        help=f"指定要并行运行的GetOrganelle进程数。\n"
             f"这决定了同时处理多少个样本。\n"
             f"默认值: {default_processes} (根据您的CPU核心数推荐)。"
    )
    # --- END NEW ---

    args = parser.parse_args()
    
    current_directory = os.getcwd()
    
    read_pairs = find_read_pairs(current_directory, args.suffix1, args.suffix2)
    
    print(f"\n🔍 找到 {len(read_pairs)} 个样本对。")
    print(f"🚀 将使用 {args.processes} 个进程并行处理，每个进程内部调用GetOrganelle时使用 {args.threads} 个线程。")
    print("-" * 60)

    # --- MODIFIED: 使用 multiprocessing.Pool 替换 for 循环 ---
    # 准备传递给 run_getorganelle 函数的参数列表
    # 每个元组对应一次函数调用，包含 (prefix, r1, r2, type, threads)
    tasks = [
        (prefix, r1, r2, args.type, args.threads) 
        for prefix, r1, r2 in read_pairs
    ]
    
    # 创建一个进程池，'with' 语句可以确保进程池在使用后被正确关闭
    with Pool(processes=args.processes) as pool:
        # pool.starmap 会将 tasks 列表中的每个元组解包作为参数传递给 run_getorganelle 函数
        # 它会阻塞直到所有任务完成
        pool.starmap(run_getorganelle, tasks)
    # --- END MODIFICATION ---
        
    print("\n🎉 所有并行任务处理完毕！")

if __name__ == "__main__":
    main()