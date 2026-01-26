import os
from Bio import SeqIO

# 根據您提供的圖表，內置葉綠體基因分類字典 (保持原始大小寫以供閱讀)
GENE_CLASSIFICATION = {
    "Self-replication": {
        "Ribosomal RNA genes": ["rrn23", "rrn16", "rrn5", "rrn4.5"],
        "Transfer RNA genes": [
            "trnA-UGC", "trnC-GCA", "trnD-GUC", "trnE-UUC", "trnF-GAA", "trnFm-CAU", "trnG-UCC",
            "trnG-GCC", "trnH-GUG", "trnI-CAU", "trnI-GAU", "trnK-UUU", "trnL-CAA", "trnL-UAA",
            "trnL-UAG", "trnM-CAU", "trnN-GUU", "trnP-UGG", "trnQ-UUG", "trnR-ACG", "trnR-UCU",
            "trnS-GCU", "trnS-GGA", "trnS-UGA", "trnT-GGU", "trnT-UGU", "trnV-GAC", "trnV-UAC",
            "trnW-CCA", "trnY-GUA"
        ],
        "Small subunit of ribosome": [
            "rps11", "rps12", "rps14", "rps15", "rps16", "rps18", "rps19", "rps2", "rps3",
            "rps4", "rps7", "rps8"
        ],
        "Large subunit of ribosome": ["rpl14", "rpl16", "rpl2", "rpl20", "rpl22", "rpl23", "rpl32", "rpl33", "rpl36"],
        "DNA-dependent RNA polymerase": ["rpoA", "rpoB", "rpoC1", "rpoC2"],
        "Translational initiation factor": ["infA"]
    },
    "Genes for photosynthesis": {
        "Subunits of photosystem I": ["psaA", "psaB", "psaC", "psaI", "psaJ", "ycf4"],
        "Subunits of photosystem II": [
            "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK",
            "psbL", "psbM", "psbT", "psbZ", "psbN"
        ],
        "Subunits of cytochrome": ["petA", "petB", "petD", "petG", "petL", "petN"],
        "Subunits of ATP synthase": ["atpA", "atpB", "atpE", "atpF", "atpH", "atpI"],
        "Large subunit of Rubisco": ["rbcL"],
        "Subunits of NADH dehydrogenase": [
            "ndhA", "ndhB", "ndhC", "ndhD", "ndhE", "ndhF", "ndhG", "ndhH", "ndhI",
            "ndhJ", "ndhK"
        ]
    },
    "Other genes": {
        "Maturase": ["matK"],
        "Envelope membrane protein": ["cemA"],
        "Subunit of acetyl-CoA": ["accD"],
        "C-type cytochrome synthesis gene": ["ccsA"],
        "Protease": ["clpP"],
        "Component of TIC complex": ["ycf1"],
        "Photosystem I assembly factor": ["ycf3"],
        "Essential unknown function": ["ycf2"],
    }
}

def analyze_genbank_files_advanced():
    """
    主函數：【高級版】
    1. 不區分基因名大小寫。
    2. 報告在參考列表之外發現的新基因。
    3. 將總結報告寫入文本文件。
    """
    genbank_dir = "input_gb_files"
    output_report_file = "gene_content_report.tsv"
    summary_report_file = "gene_summary_report.txt" # 新增的總結報告文件名

    # --- 1. 準備工作：創建小寫版本的參考列表 ---
    GENE_CLASSIFICATION_LOWER = {
        cat: {group: [gene.lower() for gene in genes] for group, genes in groups.items()}
        for cat, groups in GENE_CLASSIFICATION.items()
    }
    all_reference_genes_lower = set()
    for category, groups in GENE_CLASSIFICATION_LOWER.items():
        for group_name, genes in groups.items():
            for gene in genes:
                all_reference_genes_lower.add(gene)

    # --- 2. 檢查輸入文件夾 ---
    if not os.path.isdir(genbank_dir):
        print(f"錯誤：找不到輸入文件夾 '{genbank_dir}'。")
        return
    gb_files = [f for f in os.listdir(genbank_dir) if f.endswith((".gb", ".gbk", ".genbank"))]
    if not gb_files:
        print(f"錯誤：在 '{genbank_dir}' 文件夾中沒有找到任何 GenBank 文件。")
        return
    print(f"找到 {len(gb_files)} 個 GenBank 文件，開始進行分析...")

    # --- 3. 從每個 GenBank 文件中提取基因 (統一轉為小寫) ---
    presence_data = {}
    for filename in sorted(gb_files):
        filepath = os.path.join(genbank_dir, filename)
        genes_in_file = set()
        try:
            for record in SeqIO.parse(filepath, "genbank"):
                for feature in record.features:
                    if feature.type == "gene" and "gene" in feature.qualifiers:
                        # **核心修改**: 統一轉為小寫
                        gene_name = feature.qualifiers["gene"][0].strip().lower()
                        genes_in_file.add(gene_name)
            presence_data[filename] = genes_in_file
            print(f" -> 已處理: {filename} (找到 {len(genes_in_file)} 個基因)")
        except Exception as e:
            print(f" -> 處理文件 {filename} 時發生錯誤: {e}")

    # --- 4. 生成 TSV 報告文件 ---
    print(f"\n正在生成 TSV 報告文件: '{output_report_file}'...")
    with open(output_report_file, 'w', encoding='utf-8') as f_out:
        header = ["Category", "Gene group", "Gene name"] + sorted(presence_data.keys())
        f_out.write("\t".join(header) + "\n")

        for category, groups in GENE_CLASSIFICATION.items():
            for group_name, genes in groups.items():
                for gene in sorted(genes):
                    row_data = [category, group_name, gene]
                    gene_lower = gene.lower() # 使用小寫版本進行比較
                    for filename in sorted(presence_data.keys()):
                        if gene_lower in presence_data.get(filename, set()):
                            row_data.append("+")
                        else:
                            row_data.append("-")
                    f_out.write("\t".join(row_data) + "\n")
    print(" -> TSV 報告文件已成功生成！")

    # --- 5. 生成詳細的書面總結報告 (TXT文件) ---
    print(f"正在生成書面總結報告: '{summary_report_file}'...")
    with open(summary_report_file, 'w', encoding='utf-8') as f_summary:
        f_summary.write("Gene Content Analysis Summary\n")
        f_summary.write("="*30 + "\n\n")

        # 5.1 總結缺失的基因
        f_summary.write("--- Missing Genes (Present in Reference, Absent in Sample) ---\n")
        missing_summary = {}
        for gene in sorted(list(all_reference_genes_lower)):
            missing_in_files = [os.path.splitext(fn)[0] for fn, found in presence_data.items() if gene not in found]
            if missing_in_files:
                missing_summary[gene] = missing_in_files
        
        if not missing_summary:
            f_summary.write("No missing genes found across all samples based on the reference list.\n")
        else:
            for gene, files in missing_summary.items():
                if len(files) == len(gb_files):
                    f_summary.write(f"Gene: {gene:<10}\tMISSING in ALL {len(files)} samples.\n")
                else:
                    f_summary.write(f"Gene: {gene:<10}\tMissing in {len(files)} sample(s): {', '.join(files)}\n")
        
        # 5.2 總結額外發現的基因
        f_summary.write("\n\n--- Extra Genes Found (Present in Sample, Absent in Reference) ---\n")
        all_found_genes = set.union(*presence_data.values())
        extra_genes = sorted(list(all_found_genes - all_reference_genes_lower))

        if not extra_genes:
            f_summary.write("No extra genes were found that are not on the reference list.\n")
        else:
            for gene in extra_genes:
                found_in_files = [os.path.splitext(fn)[0] for fn, found in presence_data.items() if gene in found]
                f_summary.write(f"Gene: {gene:<15}\tFound in {len(found_in_files)} sample(s): {', '.join(found_in_files)}\n")
    
    print(f" -> 書面總結報告已成功生成！")
    print("\n所有分析已完成。")

if __name__ == "__main__":
    analyze_genbank_files_advanced()