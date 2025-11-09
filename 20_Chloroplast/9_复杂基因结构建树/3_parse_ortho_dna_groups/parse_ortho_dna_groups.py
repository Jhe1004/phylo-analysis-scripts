import os
import glob
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# --- 脚本配置 ---

# 1. 过滤器：是否只保留单拷贝直系同源基因 (推荐用于建树)
SINGLE_COPY_ONLY = True

# 2. 过滤器：一个OG至少需要包含多少个物种才会被提取
MIN_SPECIES_PER_OG = 30

# 3. 查找TSV子集文件的通配符
TSV_GLOB_PATTERN = "proteinortho_*.tsv"

# --- 辅助函数 ---

def load_all_sequences(species_files):
    """
    从所有物种的FASTA文件中加载所有序列片段到内存中。
    返回一个字典：{ 'seq_id': SeqRecord }
    
    注意：此函数需要原始的、未分割的FASTA文件 (例如 'MW542990_Clematis_potaninii.fasta')
    与本脚本在同一目录中。
    """
    all_seq_dict = {}
    
    # 构建一个所有需要的 .fasta 文件的集合
    required_files = set()
    for f in species_files:
        if f.endswith(".fasta"):
            required_files.add(f)
        elif f.endswith(".gb") or f.endswith(".gbk"):
             # 如果TSV头部是 .gb 文件名, 我们需要 .fasta 文件
             base_name = os.path.splitext(f)[0]
             required_files.add(f"{base_name}.fasta")

    files_loaded = 0
    print(f"  正在从 {len(required_files)} 个源FASTA文件中加载所有DNA片段...")
    
    for fasta_file in required_files:
        if not os.path.exists(fasta_file):
            print(f"    警告: 未在当前目录找到源文件 '{fasta_file}'。将跳过此文件。")
            continue
            
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in all_seq_dict:
                    print(f"    警告: 发现重复的序列ID: {record.id}。将覆盖。")
                all_seq_dict[record.id] = record
            files_loaded += 1
        except Exception as e:
            print(f"    错误: 读取 {fasta_file} 时出错: {e}")
            
    print(f"  成功从 {files_loaded} 个文件中加载了 {len(all_seq_dict)} 个序列片段。\n")
    return all_seq_dict

def process_one_tsv(tsv_file, output_dir, all_sequences_db):
    """
    处理单个TSV文件，并将OG FASTA文件写入其指定的输出目录。
    """
    
    # 1. 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    print(f"  开始处理TSV文件。将把OG文件写入 '{output_dir}' 文件夹...")

    og_counter = 0
    og_written = 0
    
    # 2. 从TSV表头解析物种文件列表
    species_files_in_header = []
    try:
        with open(tsv_file, 'r') as f:
            for line in f:
                if line.startswith("# Species"):
                    header_cols = line.strip().split('\t')
                    species_files_in_header = header_cols[3:]
                    break
    except Exception as e:
        print(f"  读取TSV文件表头时出错: {e}")
        return # 跳过此文件

    if not species_files_in_header:
        print(f"  错误: 无法从 '{tsv_file}' 中解析出物种文件列表。")
        return # 跳过此文件

    # 3. 遍历TSV文件的每一行 (每一个OG)
    with open(tsv_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            cols = line.strip().split('\t')
            
            try:
                num_species = int(cols[0])
                num_genes = int(cols[1]) # "genes" = "fragments"
            except (ValueError, IndexError):
                continue

            # --- 应用过滤器 ---
            if num_species < MIN_SPECIES_PER_OG:
                continue
                
            if SINGLE_COPY_ONLY and num_species != num_genes:
                continue
            
            og_counter += 1
            og_name = f"OG_{og_counter:05d}"
            og_records = [] 

            # 遍历每个物种的列
            for i, species_file in enumerate(species_files_in_header):
                # 从文件名获取物种ID (例如: MW542990_Clematis_potaninii)
                species_id = os.path.splitext(species_file)[0]
                
                cell_data = cols[i + 3]
                
                if cell_data == '*':
                    continue
                
                seq_ids_in_cell = cell_data.split(',')
                
                if SINGLE_COPY_ONLY and len(seq_ids_in_cell) > 1:
                     continue

                for j, seq_id in enumerate(seq_ids_in_cell):
                    if seq_id not in all_sequences_db:
                        print(f"    警告: 在OG {og_name} 中找不到序列ID '{seq_id}'。")
                        continue
                    
                    original_record = all_sequences_db[seq_id]
                    
                    new_id = species_id
                    if len(seq_ids_in_cell) > 1:
                        new_id = f"{species_id}_p{j+1}"
                    
                    new_record = SeqRecord(
                        original_record.seq,
                        id=new_id,
                        description=f"original_id={seq_id}"
                    )
                    og_records.append(new_record)

            # 5. 将这个OG的FASTA文件写入磁盘
            if og_records:
                output_og_file = os.path.join(output_dir, f"{og_name}.fasta")
                SeqIO.write(og_records, output_og_file, "fasta")
                og_written += 1

    print(f"  处理完成 '{tsv_file}'。")
    print(f"  总共解析了 {og_counter} 个符合条件的OG。")
    print(f"  成功将 {og_written} 个FASTA文件写入 '{output_dir}' 文件夹。\n")


# --- 运行脚本 ---
if __name__ == "__main__":
    
    # 检查 BioPython 是否已安装
    try:
        from Bio import SeqIO
    except ImportError:
        print("错误: 未找到 BioPython 库。")
        print("请先安装 BioPython: pip install biopython")
        sys.exit(1)
        
    # --- 批量处理逻辑 ---
    
    print("开始批量解析 Proteinortho TSV 子集...")
    
    tsv_files = glob.glob(TSV_GLOB_PATTERN)
    
    if not tsv_files:
        print(f"错误: 未找到与 '{TSV_GLOB_PATTERN}' 匹配的TSV文件。")
        print("请确保 'proteinortho_cds.tsv' 等文件在当前目录中。")
        sys.exit(1)
        
    print(f"找到 {len(tsv_files)} 个TSV子集文件: {', '.join(tsv_files)}\n")
    
    # 1. 首先，我们需要从 *任意一个* TSV 文件中读取表头，
    #    来获取所有物种的 .fasta 文件列表，并加载它们（只需加载一次）
    
    temp_header_file = tsv_files[0]
    species_files_in_header = []
    try:
        with open(temp_header_file, 'r') as f:
            for line in f:
                if line.startswith("# Species"):
                    header_cols = line.strip().split('\t')
                    species_files_in_header = header_cols[3:]
                    break
    except Exception as e:
        print(f"读取 {temp_header_file} 表头时出错: {e}")
        sys.exit(1)

    if not species_files_in_header:
        print(f"错误: 无法从 {temp_header_file} 解析物种列表。")
        sys.exit(1)
        
    # 2. 一次性加载所有需要的序列到内存中
    all_sequences_db = load_all_sequences(species_files_in_header)
    if not all_sequences_db:
        print("错误: 未能从FASTA文件中加载任何序列。脚本终止。")
        sys.exit(1)

    # 3. 循环处理每一个TSV文件
    for tsv_file in tsv_files:
        print(f"--- 开始处理文件: {tsv_file} ---")
        
        # 自动生成输出文件夹名称
        # 例如: "proteinortho_cds.tsv" -> "cds_groups"
        base_name = os.path.splitext(tsv_file)[0]
        dir_suffix = base_name.replace("proteinortho_", "")
        output_dir = f"{dir_suffix}_groups"
        
        # 调用处理函数
        process_one_tsv(tsv_file, output_dir, all_sequences_db)
        
    print("--- 所有批量任务已完成 ---")