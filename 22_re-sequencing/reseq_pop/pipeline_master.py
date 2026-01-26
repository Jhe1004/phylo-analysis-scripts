#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
重测序全流程控制脚本
功能：
1. 整合 BWA、GATK、Consensus 等所有步骤
2. 支持开关各个阶段（True/False）
3. 参数集中管理，方便修改
"""

import os
import subprocess
import sys
import time

# =============================================================================
#                               参数配置区
# =============================================================================

# --- 1. 流程控制开关 ---
STEPS = {
    "step1_index": False,       # 建立参考基因组索引
    "step2_mapping": True,      # BWA 比对与标记重复 (支持多样本并行)
    "step3_calling": True,     # GATK 变异检测
    "step4_consensus": True,   # 生成共识序列 (VCF to FASTA)
    "step5_qc": True,          # 质量控制检查
    "step6_combine": True      # 合并输出结果
}

# --- 2. 基础路径与全局参数 ---
REF_GENOME = "homo_ref.fasta"
FASTQ_DIR = "."                # Fastq 文件所在目录
BAM_DIR = "./bam_output"       # BAM 输出目录
VCF_DIR = "./vcf_output"       # VCF 输出目录
CONSENSUS_DIR = "./consensus_fasta_output" # FASTA 输出目录

# --- 3. 步骤具体参数 ---
# Step 2: Mapping 参数
S2_THREADS_PER_SAMPLE = 15      # 每个样本占用的线程数
S2_PARALLEL_SAMPLES = 5         # 同时并行的样本数
S2_PLUS_TAG = "_clean_1.fq.gz"
S2_MINUS_TAG = "_clean_2.fq.gz"

# Step 3: GATK 参数
S3_PARALLEL_SAMPLES = 78        # GATK 总并发任务数 (Scatter-Gather 分片模式)
S3_GATK_THREADS = 1             # 每个任务的线程数 (建议减小，以增加并发数)

# Step 4: Consensus 参数
S4_PARALLEL_SAMPLES = 10
S4_MIN_GQ = 20
S4_MIN_DP = 3

# =============================================================================
#                               核心逻辑区
# =============================================================================

def run_cmd(cmd):
    """运行 shell 命令并检查退出码"""
    print(f"\n[执行命令] {cmd}")
    ret = subprocess.run(cmd, shell=True)
    if ret.returncode != 0:
        print(f"\n[错误] 命令执行失败，程序退出。")
        sys.exit(1)

def main():
    start_time = time.time()
    print("="*60)
    print(f"重测序流程启动时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*60)

    # 创建必要的目录
    for d in [BAM_DIR, VCF_DIR, CONSENSUS_DIR]:
        if not os.path.exists(d):
            os.makedirs(d)

    # --- Step 1: Index Reference ---
    if STEPS["step1_index"]:
        print("\n>>> [Step 1/6] Indexing Reference...")
        run_cmd(f"python 1_index_reference.py -r {REF_GENOME}")

    # --- Step 2: BWA Mapping ---
    if STEPS["step2_mapping"]:
        print("\n>>> [Step 2/6] BWA Mapping & MarkDuplicates (Parallel Mode)...")
        cmd = (
            f"python 2_bwa_map.py "
            f"-r {REF_GENOME} -f {FASTQ_DIR} -o {BAM_DIR} "
            f"-p {S2_PLUS_TAG} -m {S2_MINUS_TAG} "
            f"-t {S2_THREADS_PER_SAMPLE} -j {S2_PARALLEL_SAMPLES} "
            f"--use-gatk-markdup"
        )
        run_cmd(cmd)

    # --- Step 3: GATK Variant Calling ---
    if STEPS["step3_calling"]:
        print("\n>>> [Step 3/6] GATK HaplotypeCaller...")
        cmd = (
            f"python 3_gatk_call.py "
            f"-r {REF_GENOME} -b {BAM_DIR} -o {VCF_DIR} "
            f"-t {S3_PARALLEL_SAMPLES} --gatk-threads {S3_GATK_THREADS}"
        )
        run_cmd(cmd)

    # --- Step 4: VCF to Consensus ---
    if STEPS["step4_consensus"]:
        print("\n>>> [Step 4/6] Generating Consensus Sequences...")
        cmd = (
            f"python 4_vcf_to_fasta.py "
            f"-r {REF_GENOME} -v {VCF_DIR} -o {CONSENSUS_DIR} "
            f"--min-gq {S4_MIN_GQ} --min-dp {S4_MIN_DP} "
            f"-t {S4_PARALLEL_SAMPLES}"
        )
        run_cmd(cmd)

    # --- Step 5: Quality Check ---
    if STEPS["step5_qc"]:
        print("\n>>> [Step 5/6] Quality Control...")
        # 注意：脚本 5 逻辑可能需要进入目录
        cmd = f"cd {CONSENSUS_DIR} && python ../5_check_non_atcg.py | tee ../quality_report.txt"
        run_cmd(cmd)

    # --- Step 6: Combine Fasta ---
    if STEPS["step6_combine"]:
        print("\n>>> [Step 6/6] Combining Sequences...")
        run_cmd(f"python 6_combine_fasta.py")

    end_time = time.time()
    duration = (end_time - start_time) / 60
    print("\n" + "="*60)
    print(f"所有已选步骤运行完成！总耗时: {duration:.2f} 分钟")
    print("="*60)

if __name__ == "__main__":
    main()
