import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# --- 1. 定义遗传密码表 (与之前相同) ---
CODON_MAP = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}
AMINO_ACID_MAP = defaultdict(list)
for codon, aa in CODON_MAP.items():
    AMINO_ACID_MAP[aa].append(codon)

# --- 2. 定义核心计算函数 (与之前相同) ---
def calculate_rscu_for_sequence(sequence):
    codon_counts = defaultdict(int)
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            codon_counts[codon] += 1
    rscu_values = {}
    for aa, codons in AMINO_ACID_MAP.items():
        if aa == '_' or len(codons) == 1: continue
        synonymous_codon_total_count = sum(codon_counts[c] for c in codons)
        if synonymous_codon_total_count > 0:
            expected_frequency = synonymous_codon_total_count / len(codons)
            for codon in codons: rscu_values[codon] = codon_counts[codon] / expected_frequency
        else:
            for codon in codons: rscu_values[codon] = 0.0
    if codon_counts['ATG'] > 0: rscu_values['ATG'] = 1.0
    if codon_counts['TGG'] > 0: rscu_values['TGG'] = 1.0
    return rscu_values, codon_counts

# --- 3. 主分析与报告生成函数 ---
def analyze_and_generate_report():
    input_directory = "output_cds_fasta"
    output_csv_file = "rscu_matrix.csv"
    output_report_file = "codon_usage_report.txt"

    if not os.path.isdir(input_directory):
        print(f"错误：输入文件夹 '{input_directory}' 不存在。请先运行第一个脚本。")
        return

    all_species_rscu = {}
    all_species_codon_counts = {}
    total_codons_per_species = {}

    print("开始分析所有物种的密码子使用情况...")
    for filename in os.listdir(input_directory):
        if not (filename.endswith(".fasta") or filename.endswith(".fa")): continue
        filepath = os.path.join(input_directory, filename)
        record = SeqIO.read(filepath, "fasta")
        species_name = record.id
        sequence = str(record.seq)
        
        print(f"  -> 正在分析: {species_name}")
        rscu, codon_counts = calculate_rscu_for_sequence(sequence)
        all_species_rscu[species_name] = rscu
        all_species_codon_counts[species_name] = codon_counts
        total_codons_per_species[species_name] = sum(codon_counts.values())

    if not all_species_rscu:
        print("错误：未能处理任何FASTA文件。")
        return

    # --- 4. 生成RSCU矩阵并保存为CSV ---
    rscu_df = pd.DataFrame(all_species_rscu).sort_index().fillna(0)
    rscu_df.to_csv(output_csv_file)
    print(f"\nRSCU数据矩阵已保存到 '{output_csv_file}'。")

    # --- 5. 计算报告所需的统计数据 ---
    print("正在计算统计数据以生成报告...")
    
    # a) 总密码子数
    avg_total_codons = int(pd.Series(total_codons_per_species).mean())

    # b) 最高和最低使用频率的氨基酸
    total_aa_counts = defaultdict(int)
    for species_counts in all_species_codon_counts.values():
        for codon, count in species_counts.items():
            if CODON_MAP[codon] != '_': # 忽略终止密码子
                total_aa_counts[CODON_MAP[codon]] += count
    
    # 转换成pandas Series方便查找最大最小值
    total_aa_series = pd.Series(total_aa_counts)
    most_used_aa = total_aa_series.idxmax()
    least_used_aa = total_aa_series.idxmin()
    # 将氨基酸单字母缩写转换为三字母
    aa_code_map = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu','G':'Gly','H':'His',
                   'I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr',
                   'W':'Trp','Y':'Tyr','V':'Val'}
    most_used_aa_long = aa_code_map.get(most_used_aa, most_used_aa)
    least_used_aa_long = aa_code_map.get(least_used_aa, least_used_aa)

    # c) RSCU > 1 的密码子数量范围
    preferred_codon_counts = (rscu_df > 1).sum(axis=0)
    min_pref_count = preferred_codon_counts.min()
    max_pref_count = preferred_codon_counts.max()

    # d) RSCU > 1 的密码子末位碱基统计
    preferred_codons = rscu_df[rscu_df.mean(axis=1) > 1].index
    au_ending_count = sum(1 for codon in preferred_codons if codon.endswith('A') or codon.endswith('T'))
    gc_ending_count = sum(1 for codon in preferred_codons if codon.endswith('G') or codon.endswith('C'))

    # e) 最偏好和最不偏好的密码子 (基于平均RSCU值)
    avg_rscu = rscu_df.mean(axis=1)
    most_preferred_codon = avg_rscu.idxmax()
    most_preferred_val = avg_rscu.max()
    most_preferred_aa = aa_code_map.get(CODON_MAP[most_preferred_codon], CODON_MAP[most_preferred_codon])

    # 找出RSCU > 0 的密码子，再在其中找平均值最小的
    least_preferred_codon = avg_rscu[avg_rscu > 0].idxmin()
    least_preferred_val = avg_rscu[avg_rscu > 0].min()
    least_preferred_aa = aa_code_map.get(CODON_MAP[least_preferred_codon], CODON_MAP[least_preferred_codon])

    # --- 6. 生成并保存报告文件 ---
    report_content = f"""
Codon Usage Analysis Report
===========================
Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

Summary of Findings:
--------------------
- The analysis was based on the concatenated protein-coding sequences from {len(all_species_rscu)} species. The average number of codons analyzed per species was approximately {avg_total_codons}.

- Overall Amino Acid Usage:
  - The most frequently used amino acid across all species was {most_used_aa_long} ({most_used_aa}).
  - The least frequently used amino acid was {least_used_aa_long} ({least_used_aa}).

- Relative Synonymous Codon Usage (RSCU) Patterns:
  - The number of codons with an RSCU value greater than 1 (indicating preferred usage) ranged from {min_pref_count} to {max_pref_count} across the analyzed species.
  - Methionine (Met) and Tryptophan (Trp) show no usage bias (RSCU = 1.0) by definition as they are encoded by a single codon.

- Base Composition of Preferred Codons:
  - Among codons with an average RSCU > 1, {au_ending_count} codons end with A or U (T), while {gc_ending_count} codons end with G or C. This suggests a bias towards A/U-ending codons.

- Most and Least Preferred Codons (based on average RSCU across all species):
  - The most preferred codon is {most_preferred_codon} (encoding {most_preferred_aa}), with the highest average RSCU value of {most_preferred_val:.3f}.
  - The least preferred codon (among those that are used) is {least_preferred_codon} (encoding {least_preferred_aa}), with the lowest average RSCU value of {least_preferred_val:.3f}.

NOTE:
This report provides statistical summaries. The conclusion about phylogenetic clustering (e.g., "clustering tree was largely grouped into three major blocks") can only be drawn by visually inspecting the heatmap generated in the next step (Script 3).
"""
    with open(output_report_file, 'w', encoding='utf-8') as f:
        f.write(report_content)

    print(f"统计报告已成功保存到 '{output_report_file}'。")
    print("\n--- 报告内容预览 ---")
    print(report_content)


if __name__ == "__main__":
    analyze_and_generate_report()