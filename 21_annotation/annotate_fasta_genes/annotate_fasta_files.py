#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
脚本名称：annotate_fasta_files_parallel_config_v2.py
功能：
    1. 使用拟南芥蛋白质编码序列创建BLAST数据库。
    2. 并行处理指定文件夹中符合特定模式的蛋白质FASTA文件。
    3. **新增**: 在运行BLAST前，清理输入蛋白质序列 (移除 '-' 和 '*')。
    4. 对每个清理后的FASTA文件，在独立的临时目录中运行BLASTP比对。
    5. 根据BLASTP结果确定最相似的拟南芥基因 (移除版本号)。
    6. 重命名原始FASTA文件为该基因名称，并保存到指定输出文件夹中。
    7. 同时找到并重命名对应的（符合指定模式的）DNA序列FASTA文件。
    8. 将重命名后的原始FASTA文件和DNA序列文件一起保存到输出文件夹中。
    9. 清理临时文件和目录。

依赖：
    - BLAST+ (makeblastdb, blastp)
    - Python 3
    - Biopython (`pip install biopython`)
    - Pandas (`pip install pandas`)

使用方法：
    1. 修改下面的 "--- Script Configuration ---" 部分。
    2. 确保 Biopython 已安装。
    3. 直接运行脚本: python annotate_fasta_files_parallel_config_v2.py
"""

import os
import subprocess
import pandas as pd
import glob
import shutil
import multiprocessing
import tempfile
import functools
import time
from Bio import SeqIO  # <--- 需要 Biopython

# ==============================================================================
# --- Script Configuration ---
# ==============================================================================
# --- Paths ---
INPUT_DIR = r"/home/hejian/work/transcriptome_pyloneny/31_annotation_SCOGs/raw"  # <--- 从错误日志更新
ARABIDOPSIS_FASTA = r"./tair_filter.fasta.transdecoder.pep" # !! 请确保此路径正确 !!
OUTPUT_DIR = r"/home/hejian/work/transcriptome_pyloneny/31_annotation_SCOGs/annotated_fasta_files" # <--- 从错误日志更新

# --- Filename Patterns ---
PEP_FILENAME_PATTERN = "*_pep_maffted.fas"
PEP_SUFFIX = "_pep_maffted.fas"
CDS_FILENAME_TEMPLATE = "{}_cds_maffted.fas"

# --- BLAST Settings ---
DB_NAME = "arabidopsis_db"
EVALUE_THRESHOLD = 1e-5
MAX_TARGET_SEQS = 1

# --- Performance Settings ---
NUM_PROCESSES = 0 # 0: auto (cores/2), 1: single, N: N processes
THREADS_PER_BLAST = 1 # 推荐为 1 当 NUM_PROCESSES > 1
# ==============================================================================
# --- End Configuration ---
# ==============================================================================


def make_blast_db(arabidopsis_fasta, db_name):
    """使用拟南芥FASTA文件创建BLAST数据库。"""
    db_name_abs = os.path.abspath(db_name)
    arabidopsis_fasta_abs = os.path.abspath(arabidopsis_fasta)
    db_files = [f"{db_name_abs}.{ext}" for ext in ['phr', 'pin', 'psq']]
    if all(os.path.exists(f) for f in db_files):
        print(f"BLAST数据库 '{db_name_abs}' 已存在，跳过创建。")
        return True
    cmd = ["makeblastdb", "-in", arabidopsis_fasta_abs, "-dbtype", "prot", "-out", db_name_abs]
    print(f"正在使用 '{arabidopsis_fasta_abs}' 创建BLAST数据库 '{db_name_abs}'...")
    try:
        process = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        print(f"BLAST数据库 '{db_name_abs}' 创建完成。")
        return True
    except FileNotFoundError:
        print("错误：找不到 'makeblastdb' 命令。请确保 BLAST+ 已安装并添加到系统 PATH。")
        return False
    except subprocess.CalledProcessError as e:
        print(f"错误：创建BLAST数据库失败。")
        print(f"命令: {' '.join(e.cmd)}")
        print(f"返回码: {e.returncode}")
        print(f"错误输出:\n{e.stderr}")
        return False
    except Exception as e:
        print(f"创建BLAST数据库时发生未知错误: {e}")
        return False

# run_blastp_in_tempdir 函数现在接收清理后的查询文件路径
def run_blastp_in_tempdir(cleaned_query_fasta_abs, db_name_abs, temp_blast_output_abs, evalue=1e-5, max_target_seqs=1, num_threads=1):
    """在当前工作目录（应该是临时目录）运行BLASTP。"""
    cmd = [
        "blastp", "-query", cleaned_query_fasta_abs, # <--- 使用清理后的文件
        "-db", db_name_abs, "-out", temp_blast_output_abs,
        "-evalue", str(evalue),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-max_target_seqs", str(max_target_seqs), "-num_threads", str(num_threads)
    ]
    try:
        process = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        return True
    except FileNotFoundError:
        print(f"错误 [进程 {os.getpid()}]: 找不到 'blastp' 命令。")
        return False
    except subprocess.CalledProcessError as e:
        # 打印更详细的错误信息，帮助调试
        print(f"错误 [进程 {os.getpid()}]: BLASTP 比对失败 for query '{os.path.basename(cleaned_query_fasta_abs)}' (源文件可能不同)。")
        print(f"  命令: {' '.join(e.cmd)}")
        print(f"  返回码: {e.returncode}")
        try:
            stderr_decoded = e.stderr
        except Exception:
             stderr_decoded = repr(e.stderr)
        print(f"  错误输出:\n{stderr_decoded}")
        return False
    except Exception as e:
        print(f"错误 [进程 {os.getpid()}]: 运行 BLASTP 时发生未知错误 for query '{os.path.basename(cleaned_query_fasta_abs)}': {e}")
        return False

def parse_blast_output(blast_output):
    """解析BLAST输出，返回一个DataFrame。"""
    columns = ['query_id', 'subject_id', 'pident', 'length', 'mismatch',
               'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
        return None
    try:
        blast_df = pd.read_csv(blast_output, sep='\t', names=columns)
        if blast_df.empty:
            return None
        return blast_df
    except pd.errors.EmptyDataError:
         return None
    except Exception as e:
        print(f"错误 [进程 {os.getpid()}]: 无法读取或解析BLAST输出文件 '{blast_output}'。\n{e}")
        return None


def determine_best_gene(blast_df):
    """根据BLAST结果确定最相似的拟南芥基因。"""
    if blast_df is None or blast_df.empty:
        return None
    try:
        if 'subject_id' not in blast_df.columns or 'bitscore' not in blast_df.columns:
             print(f"错误 [进程 {os.getpid()}]: BLAST DataFrame 缺少 'subject_id' 或 'bitscore' 列。")
             return None
        gene_scores = blast_df.groupby('subject_id')['bitscore'].sum().reset_index()
        best_gene_row = gene_scores.sort_values(by='bitscore', ascending=False).iloc[0]
        best_gene = best_gene_row['subject_id']
        # 移除可能的版本号，例如 AT1G01010.1 -> AT1G01010
        if '.' in best_gene:
            best_gene = best_gene.split('.')[0]
        return best_gene
    except IndexError:
        print(f"警告 [进程 {os.getpid()}]: 无法确定最佳基因（可能分组后为空）。")
        return None
    except Exception as e:
        print(f"错误 [进程 {os.getpid()}]: 在确定最佳基因时发生错误: {e}")
        return None


def find_corresponding_cds_file(pep_fasta_path, input_dir_abs, pep_suffix, cds_template):
    """根据蛋白质FASTA文件路径、其后缀和CDS模板查找对应的CDS FASTA文件路径。"""
    base_name = os.path.basename(pep_fasta_path)
    if not base_name.endswith(pep_suffix):
        return None
    stem = base_name[:-len(pep_suffix)]
    if not stem:
         return None
    expected_cds_filename = cds_template.format(stem)
    cds_search_path = os.path.join(input_dir_abs, expected_cds_filename)
    if os.path.isfile(cds_search_path):
        return cds_search_path
    else:
        return None

def clean_fasta_for_blast(original_fasta_path, cleaned_fasta_path):
    """
    读取原始FASTA，移除'-'和'*'，写入新的清理过的FASTA文件。
    返回 True 如果成功, False 如果失败。
    """
    try:
        cleaned_records = []
        for record in SeqIO.parse(original_fasta_path, "fasta"):
            # 获取序列字符串并清理
            seq_string = str(record.seq)
            cleaned_seq_string = seq_string.replace("-", "").replace("*", "")
            # 更新记录的序列
            record.seq = type(record.seq)(cleaned_seq_string) # 保持原始序列类型
            cleaned_records.append(record)

        # 写入清理后的记录到新文件
        with open(cleaned_fasta_path, "w") as outfile:
            SeqIO.write(cleaned_records, outfile, "fasta")
        return True
    except Exception as e:
        print(f"错误 [进程 {os.getpid()}]: 清理 FASTA 文件 '{os.path.basename(original_fasta_path)}' 时失败: {e}")
        return False


def process_single_fasta(original_pep_fasta_path, # 变化参数
                         # --- 以下为固定参数 ---
                         input_dir_abs, output_dir_abs, db_name_abs,
                         pep_suffix, cds_template,
                         evalue, max_target_seqs, num_threads_per_blast):
    """处理单个蛋白质FASTA文件的完整流程。"""
    process_id = os.getpid()
    base_name = os.path.basename(original_pep_fasta_path)
    print(f"[进程 {process_id}] 开始处理: {base_name}")
    temp_dir = None
    best_gene_found = None
    status = "失败"
    error_message = ""

    try:
        # 1. 创建临时目录
        temp_dir_parent = os.path.dirname(output_dir_abs) # 在输出目录旁边创建临时文件夹
        temp_dir = tempfile.mkdtemp(prefix=f'blast_{process_id}_{base_name}_', dir=temp_dir_parent)

        # 2. 清理输入 FASTA 文件 (移除 - 和 *)
        cleaned_fasta_path = os.path.join(temp_dir, f"cleaned_{base_name}")
        if not clean_fasta_for_blast(original_pep_fasta_path, cleaned_fasta_path):
            error_message = "FASTA 文件清理失败"
            raise RuntimeError(error_message)

        # 3. 定义临时 BLAST 输出文件路径
        temp_blast_output = os.path.join(temp_dir, f"{base_name}.blastp.out")

        # 4. 运行 BLASTP (使用清理后的 FASTA 文件)
        blast_success = run_blastp_in_tempdir(cleaned_fasta_path, db_name_abs, temp_blast_output, evalue, max_target_seqs, num_threads_per_blast)
        if not blast_success:
            error_message = "BLASTP 运行失败 (可能输入文件仍有问题或BLAST本身错误)"
            raise RuntimeError(error_message)

        # 5. 解析 BLAST 输出
        blast_df = parse_blast_output(temp_blast_output)
        if blast_df is None:
            print(f"  [进程 {process_id}] 文件 '{base_name}' 没有 BLAST 命中或无法解析结果。")
            status = "无命中"
            # 注意: 即使无命中，也返回原始文件名等信息
            return base_name, status, None, "无 BLAST 命中或解析错误"

        # 6. 确定最佳匹配基因
        best_gene = determine_best_gene(blast_df)
        if not best_gene:
            print(f"  [进程 {process_id}] 无法确定文件 '{base_name}' 的最佳匹配基因。")
            status = "无最佳基因"
            return base_name, status, None, "无法确定最佳基因"

        best_gene_found = best_gene
        print(f"  [进程 {process_id}] 文件 '{base_name}' 的最佳匹配基因: '{best_gene}'")

        # 7. 准备重命名和复制原始文件 (不是清理后的文件)
        _, pep_ext = os.path.splitext(base_name)
        new_pep_file_name = f"{best_gene}_pep{pep_ext}"
        final_output_pep_path = os.path.join(output_dir_abs, new_pep_file_name)

        # 查找对应的原始CDS文件
        original_cds_path = find_corresponding_cds_file(original_pep_fasta_path, input_dir_abs, pep_suffix, cds_template)

        # 8. 复制并重命名原始蛋白质文件
        try:
            # 复制原始文件，而不是清理后的文件
            shutil.copy2(original_pep_fasta_path, final_output_pep_path)
        except Exception as e:
            error_message = f"复制蛋白质文件 '{base_name}' 失败: {e}"
            raise RuntimeError(error_message)

        # 9. 如果找到CDS文件，复制并重命名它
        if original_cds_path:
            cds_base_name = os.path.basename(original_cds_path)
            _, cds_ext = os.path.splitext(cds_base_name)
            new_cds_file_name = f"{best_gene}_cds{cds_ext}"
            final_output_cds_path = os.path.join(output_dir_abs, new_cds_file_name)
            try:
                shutil.copy2(original_cds_path, final_output_cds_path)
            except Exception as e:
                print(f"  [进程 {process_id}] 警告: 复制CDS文件 '{cds_base_name}' 失败: {e}")
                error_message += f"; CDS复制失败: {e}" # 附加错误信息
        else:
            print(f"  [进程 {process_id}] 警告: 未找到 '{base_name}' 对应的CDS文件，跳过CDS复制。")
            error_message += "; 未找到CDS文件" # 附加信息

        status = "成功"

    except Exception as e:
        status = "失败"
        if not error_message: error_message = str(e)
        print(f"错误 [进程 {process_id}]: 处理文件 '{base_name}' 时发生错误: {error_message}")

    finally:
        # 10. 清理临时目录
        if temp_dir and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                print(f"警告 [进程 {process_id}]: 无法删除临时目录 '{temp_dir}': {e}")

    print(f"[进程 {process_id}] 完成处理: {base_name} (状态: {status})")
    # 返回原始文件名、状态、最佳基因、错误信息
    return base_name, status, best_gene_found, error_message


def run_pipeline():
    """封装了主流程的函数。"""
    # --- 使用配置变量 ---
    input_dir = INPUT_DIR
    arabidopsis_fasta = ARABIDOPSIS_FASTA
    output_dir = OUTPUT_DIR
    db_name = DB_NAME
    pep_pattern = PEP_FILENAME_PATTERN
    pep_suffix = PEP_SUFFIX
    cds_template = CDS_FILENAME_TEMPLATE
    evalue = EVALUE_THRESHOLD
    max_target_seqs = MAX_TARGET_SEQS
    num_processes_config = NUM_PROCESSES
    threads_per_blast = THREADS_PER_BLAST

    # --- 参数验证和准备 ---
    print("--- 步骤 0: 检查配置和环境 ---")
    input_dir_abs = os.path.abspath(input_dir)
    arabidopsis_fasta_abs = os.path.abspath(arabidopsis_fasta)
    output_dir_abs = os.path.abspath(output_dir)
    db_name_abs = os.path.abspath(db_name)

    # (路径和模式检查基本同前一版本)
    # ... [省略部分重复的检查代码，假设它们存在且有效] ...
    if not os.path.isdir(input_dir_abs): print(f"错误：输入目录 '{input_dir_abs}' 无效。"); return
    if not os.path.isfile(arabidopsis_fasta_abs): print(f"错误：拟南芥 FASTA '{arabidopsis_fasta_abs}' 无效。"); return
    if not pep_pattern or not pep_suffix or not cds_template: print(f"错误：文件名模式配置不完整。"); return
    if '{}' not in cds_template: print(f"错误：CDS_FILENAME_TEMPLATE 必须包含 '{{}}'。"); return

    # 检查 Biopython 是否导入成功
    try:
        SeqIO
    except NameError:
        print("错误：Biopython 库未找到或导入失败。请运行 'pip install biopython'。")
        return

    # 创建输出目录 (如果需要)
    # ... [省略目录创建代码] ...
    if not os.path.exists(output_dir_abs): os.makedirs(output_dir_abs)
    elif not os.path.isdir(output_dir_abs): print(f"错误: 输出路径 '{output_dir_abs}' 不是目录。"); return
    print(f"输出目录 '{output_dir_abs}' 已准备就绪。")
    print("-" * 30)

    # --- 创建BLAST数据库 ---
    print("--- 步骤 1: 准备 BLAST 数据库 ---")
    if not make_blast_db(arabidopsis_fasta_abs, db_name_abs): return
    print("-" * 30)

    # --- 获取要处理的文件列表 ---
    print(f"--- 步骤 2: 查找输入文件 (模式: {pep_pattern}) ---")
    search_pattern_abs = os.path.join(input_dir_abs, pep_pattern)
    pep_fasta_files = glob.glob(search_pattern_abs)
    if not pep_fasta_files: print(f"警告：未找到匹配 '{pep_pattern}' 的文件。"); return
    num_files = len(pep_fasta_files)
    print(f"找到 {num_files} 个匹配的蛋白质FASTA文件。")
    print("-" * 30)

    # --- 设置并行处理 ---
    print("--- 步骤 3: 并行处理文件 ---")
    # (计算 num_workers 同前一版本)
    if num_processes_config == 0:
        cpu_count = os.cpu_count()
        num_workers = cpu_count // 2 if cpu_count and cpu_count > 1 else 1
    else:
        num_workers = num_processes_config
    num_workers = max(1, min(num_workers, num_files))

    if num_workers == 1: print("以单进程模式运行。")
    else: print(f"将使用 {num_workers} 个并行进程。")
    print(f"每个 BLASTP 任务将使用 {threads_per_blast} 个线程。")
    # (CPU 过载警告同前一版本)
    # ...
    print("-" * 30)

    # --- 使用 multiprocessing.Pool 执行任务 ---
    start_time = time.time()
    results = []
    worker_func = functools.partial(
        process_single_fasta, # 注意这里函数名没有变
        # 固定参数:
        input_dir_abs=input_dir_abs, output_dir_abs=output_dir_abs, db_name_abs=db_name_abs,
        pep_suffix=pep_suffix, cds_template=cds_template,
        evalue=evalue, max_target_seqs=max_target_seqs,
        num_threads_per_blast=threads_per_blast
    )

    if num_workers > 1:
        try:
            with multiprocessing.Pool(processes=num_workers) as pool:
                results = pool.map(worker_func, pep_fasta_files)
        except Exception as e:
             print(f"\n错误：运行进程池时发生严重错误: {e}"); return
    else:
        # 单进程模式
        for fasta_file in pep_fasta_files:
            results.append(worker_func(fasta_file))

    end_time = time.time()
    print("-" * 30)

    # --- 报告结果 ---
    print("--- 步骤 4: 处理完成，生成摘要 ---")
    # (结果统计代码同前一版本)
    successful_count = 0; no_hit_count = 0; no_gene_count = 0; failed_count = 0
    failed_files = []
    for result in results:
        if isinstance(result, tuple) and len(result) == 4:
            base_name, status, best_gene, error_msg = result
            if status == "成功": successful_count += 1
            elif status == "无命中": no_hit_count += 1
            elif status == "无最佳基因": no_gene_count += 1
            else:
                failed_count += 1
                failed_files.append((base_name, error_msg))
        else:
             failed_count += 1
             failed_files.append(("未知文件/格式错误", f"工作函数返回了意外结果: {result}"))

    print(f"总共尝试处理文件数: {num_files}") # 改为“尝试处理”
    print(f"成功注释并重命名: {successful_count}")
    print(f"无 BLAST 命中或解析错误: {no_hit_count}")
    print(f"无法确定最佳基因: {no_gene_count}")
    print(f"处理失败 (FASTA清理/BLAST/复制等错误): {failed_count}") # 描述更清晰

    if failed_files:
        print("\n失败的文件列表及原因:")
        for fname, msg in failed_files:
            print(f"  - {fname}: {msg}")

    print(f"\n总耗时: {end_time - start_time:.2f} 秒")
    print("脚本执行完毕。")


if __name__ == "__main__":
    run_pipeline()