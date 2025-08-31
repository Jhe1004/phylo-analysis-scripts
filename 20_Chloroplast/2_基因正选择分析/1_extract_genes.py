import os
import re
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import defaultdict

def clean_name(name):
    """清理物種名稱。"""
    special_chars = r" |\(|\)|\[|\]|\.|\,"
    cleaned = re.sub(special_chars, "_", name)
    cleaned = re.sub(r"__+", "_", cleaned)
    return cleaned.strip('_')

def extract_shared_genes():
    """
    主函數，從所有 GenBank 檔案中提取共享基因的 CDS 和 PEP 序列。
    【偵錯強化版】能捕捉 TranslationError 並報告問題基因。
    """
    input_directory = "input_gb_files"
    output_cds_dir = "output_cds_genes"
    output_pep_dir = "output_pep_genes"

    if not os.path.isdir(input_directory):
        print(f"錯誤：輸入文件夾 '{input_directory}' 不存在。")
        return
    for d in [output_cds_dir, output_pep_dir]:
        if not os.path.exists(d):
            os.makedirs(d)

    print("第一輪：掃描所有物種以確定共享基因...")
    gene_counts = defaultdict(int)
    all_species_genes = {}
    species_name_map = {}

    gb_files = [f for f in os.listdir(input_directory) if f.endswith(".gb") or f.endswith(".gbk")]
    num_species = len(gb_files)

    for filename in gb_files:
        original_base_name = os.path.splitext(filename)[0]
        species_name = clean_name(original_base_name)
        species_name_map[species_name] = filename
        all_species_genes[species_name] = {}
        filepath = os.path.join(input_directory, filename)
        record = SeqIO.read(filepath, "genbank")
        
        for feature in record.features:
            if feature.type == "CDS" and "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
                if gene_name in ["ycf1", "rps12"]: continue
                gene_counts[gene_name] += 1
                all_species_genes[species_name][gene_name] = feature

    shared_genes = {gene for gene, count in gene_counts.items() if count == num_species}
    print(f"掃描完成。共發現 {len(shared_genes)} 個共享基因。")
    if not shared_genes:
        return

    print("\n第二輪：提取並寫入共享基因的 CDS 和蛋白質序列...")
    successful_genes = 0
    for gene in sorted(list(shared_genes)):
        cds_records = []
        pep_records = []
        
        try:
            for species_name in sorted(species_name_map.keys()):
                original_filename = species_name_map[species_name]
                original_filepath = os.path.join(input_directory, original_filename)
                
                feature_location = all_species_genes[species_name][gene].location
                full_record = SeqIO.read(original_filepath, "genbank")
                sequence_record = feature_location.extract(full_record)

                sequence_record.id = species_name
                sequence_record.description = f"[gene={gene}]"
                cds_records.append(sequence_record)
                
                # **【核心修改】** 將翻譯操作放在 try...except 塊中
                # 這樣我們可以在出錯時知道是哪個物種和基因
                protein_seq = sequence_record.seq.translate(table=11, cds=True)
                protein_record = sequence_record[:]
                protein_record.seq = protein_seq
                pep_records.append(protein_record)

            # 如果所有物種都成功翻譯，則寫入檔案
            SeqIO.write(cds_records, os.path.join(output_cds_dir, f"{gene}.fasta"), "fasta")
            SeqIO.write(pep_records, os.path.join(output_pep_dir, f"{gene}.fasta"), "fasta")
            successful_genes += 1

        except CodonTable.TranslationError as e:
            # **【核心修改】** 捕捉錯誤並打印詳細報告
            print("-" * 50)
            print(f"!!! 警告：在處理基因 '{gene}' 時遇到翻譯錯誤。")
            print(f"    -> 物種: {species_name}")
            print(f"    -> 錯誤訊息: {e}")
            print(f"    -> 該基因的前15個鹼基是: {sequence_record.seq[:15]}")
            print(f"    -> 將會跳過基因 '{gene}'，不為其生成任何輸出檔案。")
            print("-" * 50)
            # 使用 continue 跳過這個 for 迴圈的當前迭代 (即跳過這個基因)
            continue
        
    print(f"\n處理完成！共成功提取了 {successful_genes} 個基因的序列。")

if __name__ == "__main__":
    extract_shared_genes()