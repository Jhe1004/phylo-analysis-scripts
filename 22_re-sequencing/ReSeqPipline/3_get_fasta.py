# -*- coding: utf-8 -*-

# 导入必要的库
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import os  # 用于文件和路径操作
import glob  # 用于查找文件模式
import multiprocessing # <--- 新增：用于并行处理

################################################################################
# 脚本参数定义 (请在此处修改您的文件路径和参数)
################################################################################

# 1. 输入参数
# 物种A的参考基因组FASTA文件路径
REF_FASTA_PATH = "ref.fasta"
# 包含一个或多个BAM文件的文件夹路径
# !!! 重要: 文件夹内所有BAM文件必须是排序过的 (sorted) 并且已经建立了索引 (indexed) !!!
# (索引文件可以是 .bai 或 .csi 后缀)
BAM_FOLDER_PATH = "./"

# 2. 输出参数
# 生成的FASTA一致性序列的输出文件夹路径
# 如果此文件夹不存在，脚本会尝试创建它
OUTPUT_FOLDER_PATH = "./consensus_fasta_output/"

# 3. 处理参数
MIN_COVERAGE_THRESHOLD = 3  # 确定一个碱基所需的最小reads覆盖度
# MIN_BASE_QUALITY_THRESHOLD 已被移除，因为原始数据没有碱基质量值
MIN_MAPPING_QUALITY_THRESHOLD = 10  # 考虑reads时的最小比对质量 (Phred score)

# 4. 并行处理参数  <--- 新增
# 设置为1将禁用并行处理 (与原始脚本行为一致)
# 通常建议设置为 os.cpu_count() 或一个适合您机器配置的数字 (例如 4, 8)
# 如果BAM文件数量少于此值，将自动使用BAM文件数量作为进程数
NUM_PARALLEL_PROCESSES = 1 # 默认使用所有CPU核心数

################################################################################
# 辅助函数 (用于Python 3.9之前的版本)
################################################################################
def custom_removesuffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[:-len(suffix)]
    return input_string

################################################################################
# 主要功能函数 (generate_consensus_fasta 基本保持不变)
################################################################################
def generate_consensus_fasta(ref_fasta_path,
                             bam_file_path,
                             output_fasta_path,
                             sample_name_for_id,
                             min_coverage,
                             min_mapping_quality):
    """
    根据单个BAM文件，为物种B生成FASTA格式的一致性基因组序列。
    此版本禁用了基于碱基质量的过滤。
    (函数内容与原脚本一致)
    """
    try:
        bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    except ValueError as e:
        # 在worker中打印文件名会更有帮助，这里保留基础错误信息
        print(f"  错误: 打开BAM文件 '{os.path.basename(bam_file_path)}' 失败: {e}")
        print(f"  请确保BAM文件 '{os.path.basename(bam_file_path)}' 已正确排序并已建立索引 (.bai 或 .csi)。")
        return False
    
    output_records = []

    # 这条打印现在由worker函数在调用前处理，或者保留但会因并行而交错
    # print(f"  正在为样本 '{sample_name_for_id}' (源自 {os.path.basename(bam_file_path)}) 生成一致性序列 (碱基质量过滤已禁用)...")

    for ref_record in SeqIO.parse(ref_fasta_path, "fasta"):
        ref_name = ref_record.id
        ref_len = len(ref_record.seq)
        consensus_b_seq_list = ['N'] * ref_len

        try:
            for pileupcolumn in bamfile.pileup(contig=ref_name,
                                               stepper="all",
                                               min_base_quality=0, # 禁用碱基质量过滤
                                               min_mapping_quality=min_mapping_quality):
                pos = pileupcolumn.pos
                if pos >= ref_len: continue

                if pileupcolumn.n >= min_coverage:
                    bases_at_pos = []
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            try:
                                base = pileupread.alignment.query_sequence[pileupread.query_position]
                                bases_at_pos.append(base.upper())
                            except IndexError:
                                pass
                    
                    if bases_at_pos:
                        base_counts = Counter(bases_at_pos)
                        most_common_base, _ = base_counts.most_common(1)[0]
                        consensus_b_seq_list[pos] = most_common_base
        except ValueError as e: # 通常是contig不在BAM header中
             # print(f"  警告: 处理参考序列 '{ref_name}' 时出现问题 (可能在BAM文件中未找到): {e}")
             pass 

        consensus_b_sequence_str = "".join(consensus_b_seq_list)
        
        record_b_id = f"{sample_name_for_id}_{ref_name}"
        record_b_description = (
            f"Consensus sequence for sample {sample_name_for_id} (derived from {os.path.basename(bam_file_path)}), "
            f"based on mapping to species A reference {ref_name}. "
            f"Parameters: min_cov={min_coverage}, min_mq={min_mapping_quality} (base quality filtering disabled)"
        )
        
        record_b = SeqRecord(Seq(consensus_b_sequence_str),
                             id=record_b_id,
                             description=record_b_description)
        output_records.append(record_b)

    if not output_records:
        print(f"  错误: 未能为BAM文件 '{os.path.basename(bam_file_path)}' 生成任何序列记录。")
        if bamfile: bamfile.close()
        return False

    SeqIO.write(output_records, output_fasta_path, "fasta")
    if bamfile: bamfile.close()
    # print(f"  一致性FASTA文件已为 '{sample_name_for_id}' 生成: {output_fasta_path}") # 由worker统一处理成功/失败信息
    return True

################################################################################
# 新增：用于并行处理的Worker函数
################################################################################
def process_bam_file_worker(args_tuple):
    """
    处理单个BAM文件的包装函数，供多进程池调用。
    """
    bam_file_path_current, ref_fasta_path_global, output_folder_path_global, \
    min_coverage_global, min_mapping_quality_global, file_index, total_files = args_tuple

    base_name_current_bam = os.path.basename(bam_file_path_current)
    sample_name = custom_removesuffix(base_name_current_bam, ".bam")
    
    print(f"[Worker {os.getpid()}] 开始处理第 {file_index + 1}/{total_files} 个BAM文件: {base_name_current_bam}")

    bam_index_bai_path = bam_file_path_current + ".bai"
    bam_index_csi_path = bam_file_path_current + ".csi"
    if not (os.path.exists(bam_index_bai_path) or os.path.exists(bam_index_csi_path)):
        print(f"  [Worker {os.getpid()}] 警告: BAM文件 '{base_name_current_bam}' 对应的索引文件 (.bai 或 .csi) 未找到。跳过此文件。")
        return (base_name_current_bam, False, "索引文件未找到")

    output_fasta_filename = f"{sample_name}_consensus.fasta"
    output_fasta_path_current = os.path.join(output_folder_path_global, output_fasta_filename)

    print(f"  [Worker {os.getpid()}] 正在为样本 '{sample_name}' (源自 {base_name_current_bam}) 生成一致性序列...")
    
    success = generate_consensus_fasta(
        ref_fasta_path=ref_fasta_path_global,
        bam_file_path=bam_file_path_current,
        output_fasta_path=output_fasta_path_current,
        sample_name_for_id=sample_name,
        min_coverage=min_coverage_global,
        min_mapping_quality=min_mapping_quality_global
    )

    if success:
        print(f"  [Worker {os.getpid()}] 成功: 为 '{sample_name}' 生成一致性FASTA: {output_fasta_path_current}")
        return (base_name_current_bam, True, output_fasta_path_current)
    else:
        print(f"  [Worker {os.getpid()}] 失败:未能为 '{sample_name}' (源自 {base_name_current_bam}) 生成一致性序列。")
        return (base_name_current_bam, False, f"处理失败 {base_name_current_bam}")

################################################################################
# 脚本执行部分
################################################################################
if __name__ == "__main__":
    # 对于Windows平台，多进程需要此保护
    multiprocessing.freeze_support() 

    print("脚本开始执行 - 批量BAM文件一致性序列生成 (多进程并行处理, 碱基质量过滤已禁用)")
    print("-" * 70)

    # 1. 检查参考FASTA文件是否存在
    if not os.path.exists(REF_FASTA_PATH):
        print(f"错误: 参考基因组FASTA文件未找到: '{REF_FASTA_PATH}'")
        exit(1)
    print(f"参考基因组: {REF_FASTA_PATH}")

    # 2. 检查并准备BAM输入文件夹和输出文件夹
    if not os.path.isdir(BAM_FOLDER_PATH):
        print(f"错误: BAM输入文件夹未找到: '{BAM_FOLDER_PATH}'")
        exit(1)
    print(f"BAM输入文件夹: {BAM_FOLDER_PATH}")

    if not os.path.isdir(OUTPUT_FOLDER_PATH):
        print(f"信息: 输出文件夹 '{OUTPUT_FOLDER_PATH}' 不存在，将自动创建。")
        try:
            os.makedirs(OUTPUT_FOLDER_PATH, exist_ok=True)
        except OSError as e:
            print(f"错误: 创建输出文件夹 '{OUTPUT_FOLDER_PATH}' 失败: {e}")
            exit(1)
    print(f"FASTA输出文件夹: {OUTPUT_FOLDER_PATH}")
    print("-" * 70)

    # 3. 查找所有BAM文件
    bam_file_pattern = os.path.join(BAM_FOLDER_PATH, "*.bam")
    bam_files_to_process = glob.glob(bam_file_pattern)

    if not bam_files_to_process:
        print(f"错误: 在文件夹 '{BAM_FOLDER_PATH}' 中未找到任何 .bam 文件。")
        exit(1)

    print(f"在文件夹 '{BAM_FOLDER_PATH}' 中找到以下 {len(bam_files_to_process)} 个BAM文件将进行处理:")
    for bf_path in bam_files_to_process:
        print(f"  - {os.path.basename(bf_path)}")
    print("-" * 70)

    # 4. 准备并行处理任务
    tasks_to_submit = []
    for i, bam_file_path_current in enumerate(bam_files_to_process):
        tasks_to_submit.append((
            bam_file_path_current,
            REF_FASTA_PATH,
            OUTPUT_FOLDER_PATH,
            MIN_COVERAGE_THRESHOLD,
            MIN_MAPPING_QUALITY_THRESHOLD,
            i, # 当前文件索引 (用于日志)
            len(bam_files_to_process) # 文件总数 (用于日志)
        ))

    successful_jobs = 0
    failed_jobs = 0
    
    if tasks_to_submit:
        if NUM_PARALLEL_PROCESSES <= 0:
            print("警告: NUM_PARALLEL_PROCESSES 设置为无效值，将默认为 1 (串行处理)。")
            effective_num_processes = 1
        else:
            effective_num_processes = min(NUM_PARALLEL_PROCESSES, len(tasks_to_submit))
        
        # 确保至少有1个进程（如果任务存在）
        if len(tasks_to_submit) > 0 and effective_num_processes == 0:
             effective_num_processes = 1

        print(f"将使用 {effective_num_processes} 个并行进程处理 {len(tasks_to_submit)} 个BAM文件。")
        print("-" * 70)

        # 使用进程池执行任务
        with multiprocessing.Pool(processes=effective_num_processes) as pool:
            results = pool.map(process_bam_file_worker, tasks_to_submit)

        # 处理结果
        for bam_name, success, message in results:
            if success:
                successful_jobs += 1
            else:
                failed_jobs += 1
                # 可以在此处打印更详细的失败原因（message中已包含）
                # print(f"  文件 {bam_name} 处理失败: {message}") 
                # (worker内部已经打印了失败信息，这里可以选择是否重复打印摘要)
    else:
        print("没有BAM文件需要处理。")


    print("-" * 70)
    print("\n所有BAM文件处理完毕。")
    print(f"总计: {len(bam_files_to_process)} 个BAM文件")
    print(f"成功处理: {successful_jobs} 个")
    print(f"失败/跳过: {failed_jobs} 个")
    print(f"FASTA输出文件保存在: {os.path.abspath(OUTPUT_FOLDER_PATH)}")
    print("-" * 70)
    print("脚本执行结束。")
    print("\n重要提示: 由于原始测序数据不包含碱基质量值，本脚本已禁用基于碱基质量的过滤。")
    print("这意味着reads中所有位置的碱基（无论其原始测序质量如何）都会被用于一致性序列的生成。")
    print("如果原始reads中包含较多测序错误，这可能会影响最终一致性序列的准确性。")