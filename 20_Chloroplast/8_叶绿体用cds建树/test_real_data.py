# -*- coding: utf-8 -*-
"""
test_all_outputs.py (Final Robust Version)

该脚本用于批量验证 extract_orthologs.py 脚本生成的所有输出文件。
它会遍历输出目录中的每个FASTA文件，并将其中每个物种的序列
与直接从输入目录中相应GenBank文件中提取的“标准答案”序列进行比较。

[最终修复]: 提取“标准答案”的逻辑与主脚本完全镜像，
优先选择 CDS > rRNA > tRNA > gene 特征，并忽略伪基因，
同时在同类型特征中优先选择最长的版本，
以确保在面对重复或有歧义的注释时（如 ycf1, rrn* 基因）验证的准确性。
"""
import re
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    print("错误：运行本测试脚本需要 Biopython 库。请使用 'pip install biopython' 命令进行安装。")
    exit(1)

# --- 配置区 ---
INPUT_DIR = Path("./input_data")
OUTPUT_DIR = Path("./output")
GB_EXTENSIONS = ['*.gb', '*.gbk']
FASTA_EXTENSION = '*.fasta'

# --- 核心辅助函数 ---

def sanitize_sample_name(name):
    """清理样本名，这个函数必须与主脚本中的函数完全一致！"""
    name = re.sub(r'[^a-zA-Z0-9_-]', '_', name)
    name = re.sub(r'__+', '_', name)
    name = name.strip('_')
    return name

def create_filename_map(input_dir):
    """创建一个从清理后的样本名到原始GenBank文件路径的映射字典。"""
    print(f"[初始化] 正在扫描输入目录 '{input_dir}' 并创建文件名映射...")
    mapping = {}
    gb_files = []
    for ext in GB_EXTENSIONS:
        gb_files.extend(list(input_dir.glob(ext)))
    
    if not gb_files:
        print(f"  - 警告：在 '{input_dir}' 中没有找到任何 GenBank 文件。")
        return None

    for gb_file in gb_files:
        sanitized_name = sanitize_sample_name(gb_file.stem)
        if sanitized_name in mapping:
            print(f"  - 警告：清理后出现重复的文件名 '{sanitized_name}'。")
        mapping[sanitized_name] = gb_file
    
    print(f"  - 完成，找到 {len(mapping)} 个唯一的 GenBank 文件进行映射。")
    return mapping

def get_best_feature(record, gene_name, candidates=None):
    """
    [最终版逻辑] 从记录或候选列表中为给定基因名选择最佳特征。
    1. 过滤掉所有伪基因。
    2. 按类型优先级 (CDS > rRNA > tRNA > gene) 和长度 (最长优先) 排序。
    3. 返回最优选。
    """
    if candidates is None:
        candidates = [f for f in record.features if 'gene' in f.qualifiers and f.qualifiers['gene'][0] == gene_name]
    
    valid_candidates = [
        f for f in candidates 
        if 'pseudo' not in f.qualifiers and 'pseudogene' not in f.qualifiers
    ]

    if not valid_candidates:
        return None, None

    priority_order = ['CDS', 'rRNA', 'tRNA', 'gene']
    
    def sort_key(feature):
        try:
            type_priority = priority_order.index(feature.type)
        except ValueError:
            type_priority = len(priority_order)
        length = len(feature)
        return (type_priority, -length)

    valid_candidates.sort(key=sort_key)
    
    best_feature = valid_candidates[0]
    extraction_type = best_feature.type
    
    return best_feature, extraction_type

def get_reference_sequence(gb_file_path, gene_name):
    """
    [最终修复逻辑] 镜像主脚本的逻辑，提取“标准答案”。
    """
    try:
        gb_record = SeqIO.read(gb_file_path, "genbank")
    except Exception as e:
        print(f"    - [错误] 无法解析GenBank文件 '{gb_file_path.name}': {e}")
        return None, None

    best_feature, extraction_type = get_best_feature(gb_record, gene_name)

    if best_feature:
        return best_feature.extract(gb_record.seq), extraction_type
    else:
        return None, None

# --- 主验证逻辑 ---

def run_verification():
    """执行完整的批量验证流程"""
    print("--- 开始批量验证所有输出的FASTA文件 ---")

    if not INPUT_DIR.is_dir() or not OUTPUT_DIR.is_dir():
        sys.exit(f"[失败] 错误：输入 '{INPUT_DIR}' 或输出 '{OUTPUT_DIR}' 目录不存在。")

    filename_map = create_filename_map(INPUT_DIR)
    if not filename_map: return
        
    fasta_files_to_test = sorted(list(OUTPUT_DIR.glob(FASTA_EXTENSION)))
    if not fasta_files_to_test:
        sys.exit(f"\n[失败] 错误：在输出目录 '{OUTPUT_DIR}' 中没有找到任何FASTA文件。")
        
    print(f"\n[信息] 准备开始测试 {len(fasta_files_to_test)} 个基因文件...\n")

    total_sequences_checked = 0
    total_passes = 0
    total_fails = 0
    failed_items = []

    for fasta_file in fasta_files_to_test:
        gene_name = fasta_file.stem
        print(f"\n--- 正在测试基因: '{gene_name}' (文件: {fasta_file.name}) ---")
        
        try:
            sequences_in_fasta = {rec.id: rec for rec in SeqIO.parse(fasta_file, "fasta")}
            if not sequences_in_fasta:
                print("  - [警告] 该FASTA文件为空，跳过。")
                continue

            # 验证所有在FASTA文件中的物种
            for species_id, record in sequences_in_fasta.items():
                total_sequences_checked += 1
                output_sequence = record.seq
                source_gb_file = filename_map.get(species_id)

                if not source_gb_file:
                    print(f"  - [FAIL] {species_id}: 无法在输入目录中找到对应的源 GenBank 文件。")
                    total_fails += 1
                    failed_items.append(f"Gene: {gene_name}, Species: {species_id} (源文件缺失)")
                    continue
                
                reference_sequence, ex_type = get_reference_sequence(source_gb_file, gene_name)

                if not reference_sequence:
                    print(f"  - [FAIL] {species_id}: 主脚本输出了序列，但测试脚本在源文件中未找到有效特征 '{gene_name}'。")
                    total_fails += 1
                    failed_items.append(f"Gene: {gene_name}, Species: {species_id} (不应被提取)")
                    continue
                
                if str(output_sequence).upper() == str(reference_sequence).upper():
                    print(f"  - [PASS] {species_id}: 序列匹配 (提取自 {ex_type})")
                    total_passes += 1
                else:
                    print(f"  - [FAIL] {species_id}: 序列不匹配！ (应提取自 {ex_type})")
                    print(f"    - 长度差异: 输出({len(output_sequence)}) vs 应为({len(reference_sequence)})")
                    total_fails += 1
                    failed_items.append(f"Gene: {gene_name}, Species: {species_id} (内容不匹配)")
            
            # 反向检查: 哪些物种本应有这个基因但没有被提取出来
            all_species_in_input = set(filename_map.keys())
            extracted_species = set(sequences_in_fasta.keys())
            missing_species = all_species_in_input - extracted_species
            
            for species_id in missing_species:
                source_gb_file = filename_map[species_id]
                reference_sequence, _ = get_reference_sequence(source_gb_file, gene_name)
                if reference_sequence: # 如果这个物种确实有这个基因
                    total_sequences_checked +=1
                    print(f"  - [FAIL] {species_id}: 该物种的基因存在于源文件中，但未被主脚本提取。")
                    total_fails += 1
                    failed_items.append(f"Gene: {gene_name}, Species: {species_id} (漏提取)")

        except Exception as e:
            print(f"  - [错误] 处理文件 '{fasta_file.name}' 时发生意外错误: {e}")
            total_fails += 1
            failed_items.append(f"Gene: {gene_name} (整个文件处理失败)")

    print("\n" + "="*40)
    print("            批量验证完成")
    print("="*40)
    print(f"总共检查的基因文件数: {len(fasta_files_to_test)}")
    print(f"总共验证的序列条目数: {total_sequences_checked}")
    print(f"           验证通过: {total_passes}")
    print(f"           验证失败: {total_fails}")
    print("="*40)

    if total_fails > 0:
        print("\n失败项目汇总：")
        for item in failed_items:
            print(f"  - {item}")
    else:
        print("\n所有序列均验证通过，恭喜！你的脚本工作得非常出色！")


if __name__ == "__main__":
    run_verification()