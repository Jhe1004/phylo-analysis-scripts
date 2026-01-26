#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import subprocess

def get_args():
    """
    获取并解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description="一个使用单进程顺序分割双端FASTQ文件的脚本 (依赖seqkit)。",
        formatter_class=argparse.RawTextHelpFormatter  # 允许在帮助信息中使用换行符
    )

    parser.add_argument(
        "-s1", "--suffix1",
        type=str,
        required=True,
        help="必须指定。第一个 (Forward/R1) 测序文件的后缀。\n"
             "例如: '_1.clean.fq.gz' 或 '_R1.fastq'"
    )
    parser.add_argument(
        "-s2", "--suffix2",
        type=str,
        required=True,
        help="必须指定。第二个 (Reverse/R2) 测序文件的后缀。\n"
             "例如: '_2.clean.fq.gz' 或 '_R2.fastq'"
    )
    
    # 定义seqkit内部线程数参数
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=4,
        help="为每个 'seqkit split2' 任务指定的线程数 (即-j参数)。\n"
             "默认值: 4。"
    )
    
    # --- NEW: 添加 -p/--parts 参数 ---
    parser.add_argument(
        "-p", "--parts",
        type=int,
        default=2,
        help="指定 'seqkit split2' 的 -p 参数值，即将文件分割成几部分。\n"
             "默认值: 2。"
    )
    # --- END NEW ---

    # 定义输出目录前缀
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="output",
        help="指定每个分割文件输出目录的前缀。\n"
             "例如: 输入文件为 'SampleA_1.fq.gz'，前缀为 'split_out'，\n"
             "则输出目录将被创建为 'split_out_SampleA'。\n"
             "默认值: 'split_output'。"
    )

    return parser.parse_args()

def get_file_list(suffix1, suffix2):
    """
    获取当前文件夹中符合目标模式的文件。
    
    :param suffix1: 用户提供的 R1 文件名后缀
    :param suffix2: 用户提供的 R2 文件名后缀
    :return: file_list: 所有匹配的 R1 文件名列表
    """
    file_list = []
    print(f"正在当前目录中搜索以 '{suffix1}' 结尾的文件...")
    for each_file in os.listdir(os.getcwd()):
        if each_file.endswith(suffix1):
            prefix = each_file[:-len(suffix1)]
            r2_file = prefix + suffix2
            
            if os.path.exists(r2_file):
                file_list.append(each_file)
            else:
                print(f"⚠️  警告: 找到 R1 文件 '{each_file}' 但未找到对应的 R2 文件 '{r2_file}'。已跳过。")
    
    if not file_list:
        print(f"❌ 错误: 未找到匹配 '{suffix1}' 的文件。请检查你的后缀字符串和文件位置。")
        exit(1)
        
    print(f"✅ 找到 {len(file_list)} 对文件待处理。")
    return file_list

# --- MODIFIED: 函数签名增加了 'parts' 参数 ---
def run_split(file_list, seqkit_jobs, output_prefix, suffix1, suffix2, parts):
    """
    对文件列表中的每个文件顺序执行 seqkit split2 命令。

    :param file_list: 所有 R1 文件名的列表
    :param seqkit_jobs: seqkit split2 使用的线程数 (-j)
    :param output_prefix: 输出目录的前缀
    :param suffix1: R1 文件的后缀
    :param suffix2: R2 文件的后缀
    :param parts: 要将文件分割成的部分数 (-p)
    """
    total_files = len(file_list)
    print(f"\n▶️ 开始顺序处理 {total_files} 对文件...")
    
    for i, r1_file in enumerate(file_list, 1):
        print("\n" + "="*50)
        print(f"处理中 ({i}/{total_files}): {r1_file}")
        
        base_name = r1_file[:-len(suffix1)]
        r2_file = base_name + suffix2
        output_dir = f"{output_prefix}_{base_name}"

        # --- MODIFIED: 在命令中使用了 'parts' 变量 ---
        command = (
            f"seqkit split2 "
            f"-1 {r1_file} "
            f"-2 {r2_file} "
            f"-p {parts} -j {seqkit_jobs} "
            f"-O ."
        )
        # --- END MODIFICATION ---
        
        print(f"  执行命令: {command}")
        try:
            result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
            print(f"✔️  成功处理 {r1_file}。")
        except subprocess.CalledProcessError as e:
            print(f"❌ 处理 {r1_file} 时发生错误: {e}")
            print(f"  错误信息 (STDERR): {e.stderr.strip()}")
        except Exception as e:
            print(f"❌ 处理 {r1_file} 时发生未知错误: {e}")

def main():
    """
    主函数
    """
    # 1. 获取参数
    args = get_args()
    
    # 2. 获取文件列表
    file_list = get_file_list(args.suffix1, args.suffix2)
    
    if not file_list:
        print("没有文件需要处理，程序退出。")
        return
        
    # --- MODIFIED: 将 'args.parts' 传递给 run_split 函数 ---
    run_split(file_list, args.jobs, args.output_prefix, args.suffix1, args.suffix2, args.parts)

    print("\n" + "="*50)
    print("🎉 所有任务已完成。")


if __name__ == "__main__":
    main()