import glob
import sys
import os
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --- 用户配置 ---
# (配置已移至函数内部或自动检测)
FILE_PATTERN = "*_maffted.fasta"
MISSING_CHAR = "?"

# --- 核心连接函数 (与上一版相同) ---
def concatenate_alignments_in_current_dir():
    """
    主函数：查找当前文件夹中所有 *_maffted.fasta 文件，
    将它们连接成一个超级矩阵，并用 '?' 填充缺失的物种。
    输出文件名将自动设置为当前文件夹的名称。
    """
    
    current_dir_name = os.path.basename(os.getcwd())
    OUTPUT_FILE = f"{current_dir_name}.fasta"
    
    print(f"  输出文件将自动命名为: {OUTPUT_FILE}")
    
    alignment_files = sorted(glob.glob(FILE_PATTERN))
    
    if not alignment_files:
        print("  错误：在此文件夹中未找到与 '{FILE_PATTERN}' 匹配的文件。跳过。\n")
        return

    print(f"  找到 {len(alignment_files)} 个比对文件。开始构建超级矩阵...")

    concatenated_seqs = {}
    all_species_ids = set()
    alignment_lengths = {}

    # 步骤 1: 扫描
    for file in alignment_files:
        try:
            alignment = AlignIO.read(file, "fasta")
            length = alignment.get_alignment_length()
            if length == 0:
                print(f"    警告: 文件 '{file}' 为空，将跳过。")
                continue
            alignment_lengths[file] = length
            for record in alignment:
                all_species_ids.add(record.id)
        except Exception as e:
            print(f"    错误: 读取 '{file}' 时出错: {e}")
            
    if not all_species_ids or not alignment_lengths:
        print("  错误: 未能从比对文件中读取任何数据。跳过。\n")
        return

    # 步骤 2: 初始化字典
    for species_id in all_species_ids:
        concatenated_seqs[species_id] = ""

    # 步骤 3: 连接
    total_length = 0
    for file in alignment_files:
        if file not in alignment_lengths:
            continue 
            
        file_length = alignment_lengths[file]
        
        seqs_in_this_file = {record.id: str(record.seq) for record in AlignIO.read(file, "fasta")}
        missing_data_segment = MISSING_CHAR * file_length
        
        for species_id in all_species_ids:
            if species_id in seqs_in_this_file:
                concatenated_seqs[species_id] += seqs_in_this_file[species_id]
            else:
                concatenated_seqs[species_id] += missing_data_segment
        
        total_length += file_length

    # 步骤 4: 写入
    output_records = []
    for species_id in sorted(all_species_ids):
        record = SeqRecord(
            Seq(concatenated_seqs[species_id]),
            id=species_id,
            description=""
        )
        output_records.append(record)

    SeqIO.write(output_records, OUTPUT_FILE, "fasta")

    print(f"  --- {current_dir_name} 处理完成 ---")
    print(f"  成功连接 {len(alignment_lengths)} 个片段。")
    print(f"  总共 {len(all_species_ids)} 个物种。")
    print(f"  超级矩阵总长度: {total_length} bp。")
    print(f"  已保存到: {OUTPUT_FILE}\n")

# --- *** 新的批量处理 'main' 逻辑 *** ---
if __name__ == "__main__":
    
    try:
        from Bio import SeqIO, AlignIO
    except ImportError:
        print("错误: 未找到 BioPython 库。")
        print("请先安装 BioPython: pip install biopython")
        sys.exit(1)
        
    main_dir = os.getcwd()
    print(f"主目录: {main_dir}")
    print("开始批量连接 (Concatenation) 任务...")
    
    # 1. 查找所有匹配的子文件夹
    #    我们假设文件夹是 *_groups (与 MAFFT 脚本一致)
    subdir_pattern = "*_groups"
    subdirectories = [d for d in glob.glob(subdir_pattern) if os.path.isdir(d)]
    
    if not subdirectories:
        print(f"错误: 在 {main_dir} 中未找到与 '{subdir_pattern}' 匹配的子文件夹。")
        sys.exit(1)
        
    print(f"找到 {len(subdirectories)} 个子文件夹: {', '.join(subdirectories)}\n")

    # 2. 遍历每个子文件夹
    for subdir in subdirectories:
        subdir_path = os.path.join(main_dir, subdir)
        print(f"--- 正在进入: {subdir} ---")
        
        try:
            os.chdir(subdir_path) # 切换到子目录
        except Exception as e:
            print(f"  错误: 无法进入文件夹 {subdir_path}: {e}\n")
            continue
        
        # 3. 在子目录中执行连接功能
        concatenate_alignments_in_current_dir()
        
        # 4. 切换回主目录
        os.chdir(main_dir)

    print("--- 所有批量连接任务已完成 ---")