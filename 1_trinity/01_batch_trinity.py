# -*- coding: utf-8 -*-
"""
脚本名称: 01_batch_trinity.py
功能: 批量扫描双端测序文件并生成/执行 Trinity 组装命令 (兼容 Singularity)
作者: Dr. Jian He (Refactored by Gemini)
依赖: Python 3.6+, Singularity
"""

import os
import sys

# ======================= ### CONFIGURATION (用户配置区) ### =======================
# 1. 输入文件设置
# 待处理文件所在的目录 (默认当前目录 ".")
INPUT_DIRECTORY = "." 

# 测序文件的后缀名列表 (脚本会按顺序尝试匹配)
# 注意：文件必须遵循 "_1" 和 "_2" 的命名规则，例如 sample_1.fq.gz
VALID_EXTENSIONS = ['.fq.gz', '.fastq.gz', '.fq', '.fastq']

# 2. Trinity 运行设置
# Singularity 镜像文件的完整路径或名称
SINGULARITY_IMAGE_PATH = "trinityrnaseq.v2.15.2.simg"

# 分配的 CPU 线程数
CPU_THREADS = 40

# 分配的最大内存 (例如 "40G", "100G")
MAX_MEMORY = "40G"

# 3. 调试设置
# 如果为 True，脚本只会打印即将执行的命令，不会真正运行 Trinity。
# 建议首次运行时设为 True 以检查命令是否正确。
DRY_RUN = False
# =================================================================================


def find_paired_samples(directory, extensions):
    """
    扫描目录，寻找符合命名规范 (_1/_2) 的成对样本。
    """
    samples = {}
    
    if not os.path.exists(directory):
        print(f"错误: 目录不存在 - {directory}")
        sys.exit(1)

    print(f"正在扫描目录: {os.path.abspath(directory)} ...")
    files_in_dir = os.listdir(directory)
    
    for ext in extensions:
        suffix_1 = f"_1{ext}"
        suffix_2 = f"_2{ext}"
        
        for filename in files_in_dir:
            if filename.endswith(suffix_1):
                base_name = filename.replace(suffix_1, "")
                file_1 = filename
                file_2 = filename.replace(suffix_1, suffix_2)
                
                if file_2 in files_in_dir:
                    if base_name not in samples:
                        samples[base_name] = {
                            'ext': ext,
                            'left': os.path.join(directory, file_1),
                            'right': os.path.join(directory, file_2)
                        }
                        print(f"  [发现样本] {base_name} (类型: {ext})")
                    else:
                        print(f"  [跳过重复] 样本 {base_name} 已存在，忽略后缀 {ext}")
    
    return samples

def run_trinity(samples):
    """
    遍历样本字典，构建并执行 Trinity 命令。
    """
    if not samples:
        print("\n未找到任何符合条件的成对文件。请检查 CONFIGURATION 中的路径和后缀设置。")
        return

    current_dir = os.getcwd()

    print("\n" + "="*50)
    print(f"准备处理 {len(samples)} 个样本")
    print("="*50 + "\n")

    for sample_id, info in samples.items():
        # 定义输出目录 (Trinity 要求输出目录不存在，或者使用 --full_cleanup)
        output_dir = os.path.join(current_dir, f"{sample_id}_trinity")
        
        # 确定 seqType (fq 或 fa)
        seq_type = "fq"
        if info['ext'] in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']:
            seq_type = "fa"

        # 构建命令
        # 注意: 使用 singularity exec 执行镜像内的 Trinity
        cmd_parts = [
            f"singularity exec -e {SINGULARITY_IMAGE_PATH} Trinity",
            f"--seqType {seq_type}",
            f"--left \"{info['left']}\"",
            f"--right \"{info['right']}\"",
            f"--CPU {CPU_THREADS}",
            f"--max_memory {MAX_MEMORY}",
            f"--output \"{output_dir}\"",
            "--full_cleanup" # 运行后清理临时文件
        ]
        
        full_command = " ".join(cmd_parts)
        
        print(f">>> 正在处理样本: {sample_id}")
        print(f"指令预览: {full_command}\n")

        if not DRY_RUN:
            exit_code = os.system(full_command)
            if exit_code == 0:
                print(f"--- 样本 {sample_id} 处理完成 ---\n")
            else:
                print(f"!!! 样本 {sample_id} 处理出错 (退出码: {exit_code}) !!!\n")
        else:
            print("--- [Dry Run] 仅打印命令，不执行 ---\n")

def main():
    samples = find_paired_samples(INPUT_DIRECTORY, VALID_EXTENSIONS)
    run_trinity(samples)

if __name__ == "__main__":
    main()