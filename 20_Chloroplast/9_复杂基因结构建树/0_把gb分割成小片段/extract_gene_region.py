import glob
import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# --- 用户配置 ---

# 1. 定义您想要作为单独片段提取的特征类型
#    (例如: 'CDS', 'rRNA', 'tRNA', 'intron', 'exon' 等)
FEATURES_TO_EXTRACT = ['CDS', 'rRNA', 'tRNA']

# --- 脚本正文 ---

def split_genome_to_granular_segments():
    """
    主函数：查找所有.gb文件, 将每个文件按特征(CDS, rRNA, tRNA)
    和间隔区(IGS)拆分。
    *** 新功能 ***:
    本脚本会正确处理 'join' 特征 (如含内含子的基因),
    将每个外显子(exon)作为单独片段提取,
    并将内含子(intron)作为单独的 'IGS' 片段提取。
    """
    
    gb_files = glob.glob("*.gb") + glob.glob("*.gbk")
    
    if not gb_files:
        print(f"错误：在当前文件夹中未找到任何 .gb 或 .gbk 文件。")
        sys.exit(1)

    print(f"总共找到 {len(gb_files)} 个 GenBank 文件。开始处理...")
    print(f"将提取的特征类型: {', '.join(FEATURES_TO_EXTRACT)}\n")

    total_files_processed = 0

    for gb_file in gb_files:
        base_name = os.path.splitext(os.path.basename(gb_file))[0]
        output_filename = f"{base_name}.fasta"
        
        try:
            record = SeqIO.read(gb_file, "genbank")
            print(f"--- 正在处理: {gb_file} (总长度: {len(record.seq)} bp) ---")

            fasta_segments = []
            current_position = 0
            total_length = len(record.seq)

            # 1. 筛选出我们感兴趣的特征
            features_of_interest = []
            for f in record.features:
                if f.type in FEATURES_TO_EXTRACT:
                    features_of_interest.append(f)
            
            # 2. *** 新逻辑：解构 'join' 特征 ***
            #    我们将把 'join' 特征 (如 rpoC1) 拆分成它们
            #    各自独立的部分 (外显子), 并将它们视为独立的特征。
            all_parts = []
            for feature in features_of_interest:
                gene_name = feature.qualifiers.get("gene", [feature.type])[0]
                feature_type = feature.type
                
                # Biopython 中, join/complement 特征的 location.parts > 1
                if len(feature.location.parts) > 1:
                    # 这是 'join' 特征 (例如 rpoC1 )
                    for i, part in enumerate(feature.location.parts):
                        # 为每个外显子创建唯一的名称
                        part_name = f"{gene_name}_exon{i+1}"
                        # 存储 (起始, 结束, 名称, 类型, SeqFeature.Location 对象)
                        all_parts.append((int(part.start), int(part.end), part_name, feature_type, part))
                else:
                    # 这是简单的连续特征 (例如 rps4 [cite: 13, 560])
                    part = feature.location
                    all_parts.append((int(part.start), int(part.end), gene_name, feature_type, part))

            # 3. 按起始位置对所有 *部分* (外显子) 进行排序
            #    这是确保正确提取间隔区和内含子的关键
            sorted_parts = sorted(all_parts, key=lambda p: p[0])

            # 4. 遍历排序后的 *部分* (外显子)，提取“间隔区”和“特征区”
            for (p_start, p_end, p_name, p_type, p_location) in sorted_parts:

                # --- a. 检查并提取间隔区 (IGS) / 内含子 ---
                if p_start > current_position:
                    # 如果特征的起始 > 游标，说明中间有间隔
                    gap_seq = record.seq[current_position:p_start]
                    
                    # 创建FASTA头部
                    gap_id = f"{base_name}_IGS_{current_position + 1}-{p_start}"
                    
                    gap_record = SeqRecord(gap_seq, id=gap_id, description="")
                    fasta_segments.append(gap_record)
                    print(f"  提取 IGS/Intron: {current_position + 1}..{p_start}")

                # --- b. 提取特征部分 (Exon) ---
                # p_location.extract() 会自动处理正/反向链
                part_seq = p_location.extract(record.seq)
                
                # 为特征创建一个描述性的ID
                feature_id = f"{base_name}_{p_name}_{p_start + 1}-{p_end}"

                feature_record = SeqRecord(part_seq, id=feature_id, description="")
                fasta_segments.append(feature_record)
                print(f"  提取 {p_type}: {p_name} ({p_start + 1}..{p_end})")

                # --- c. 更新游标位置 ---
                # 游标总是更新到当前已处理片段的末尾
                current_position = max(current_position, p_end)

            # 5. 提取最后一个间隔区 (从最后一个特征到基因组末端)
            if current_position < total_length:
                last_igs_seq = record.seq[current_position:total_length]
                last_igs_id = f"{base_name}_IGS_{current_position + 1}-{total_length}"
                last_igs_record = SeqRecord(last_igs_seq, id=last_igs_id, description="")
                fasta_segments.append(last_igs_record)
                print(f"  提取 末尾IGS: {current_position + 1}..{total_length}")

            # 6. 将此文件的所有片段写入同名的FASTA文件
            SeqIO.write(fasta_segments, output_filename, "fasta")
            print(f"  成功: 已保存 {len(fasta_segments)} 个片段到 {output_filename}\n")
            total_files_processed += 1

        except Exception as e:
            print(f"    *** 错误: 处理 {gb_file} 时遇到问题: {e} ***\n")

    print(f"处理完成！总共成功处理了 {total_files_processed} 个文件。")


# --- 运行脚本 ---
if __name__ == "__main__":
    
    # 检查 BioPython 是否已安装
    try:
        from Bio import SeqIO
    except ImportError:
        print("错误: 未找到 BioPython 库。")
        print("请先安装 BioPython: pip install biopython")
        sys.exit(1)
        
    split_genome_to_granular_segments()