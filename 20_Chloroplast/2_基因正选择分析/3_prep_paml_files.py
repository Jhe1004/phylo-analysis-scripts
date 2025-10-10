import os
import shutil
from Bio import AlignIO
from Bio import Phylo
from io import StringIO

def write_sequential_phylip(alignment, file_path):
    """
    【优化版】手動將比對物件寫入 PAML 偏好的順序式 (sequential) PHYLIP 格式。
    此版本會自動檢測最長的物種名，並以此為標準進行對齊，以適應原始物種名。
    """
    num_species = len(alignment)
    seq_length = alignment.get_alignment_length()

    # 找出所有物種ID中最長的長度
    max_id_len = 0
    for record in alignment:
        if len(record.id) > max_id_len:
            max_id_len = len(record.id)

    with open(file_path, 'w', newline='\n') as f_out:
        # 寫入 PHYLIP 表頭
        f_out.write(f" {num_species} {seq_length}\n")
        # 寫入每個物種的序列
        for record in alignment:
            # 根據最長的名字長度進行填充，確保後面的序列能夠對齊
            padded_name = record.id.ljust(max_id_len) + "  "
            f_out.write(padded_name + str(record.seq) + "\n")

def write_paml_newick(tree, file_path):
    """
    手動將樹物件寫入 PAML 專用的 Newick 格式，包含物種數量表頭。
    (此函數無需修改)
    """
    num_species = tree.count_terminals()
    handle = StringIO()
    Phylo.write(tree, handle, "newick")
    newick_string = handle.getvalue().strip()

    with open(file_path, 'w', newline='\n') as f_out:
        f_out.write(f" {num_species} 1\n")
        f_out.write(newick_string + "\n")


def prepare_paml_inputs_final_final():
    """
    主函數，準備 PAML/codeml 所需的所有輸入檔案。
    【修改版】此版本直接使用原始物種名，不再進行縮寫替換。
    """
    aligned_cds_dir = "aligned_cds"
    paml_input_dir = "paml_input"
    tree_file = "RAxML_bipartitions.result.newick" # 這是您原始的、不帶表頭的標準 Newick 檔案

    if not os.path.isdir(aligned_cds_dir) or not os.path.exists(tree_file):
        print(f"錯誤：找不到輸入文件夾 '{aligned_cds_dir}' 或樹檔案 '{tree_file}'。")
        return

    if os.path.exists(paml_input_dir):
        shutil.rmtree(paml_input_dir)
    os.makedirs(paml_input_dir)

    print("讀取原始樹檔案...")
    try:
        # 直接讀取樹對象，不再進行名稱替換
        tree = Phylo.read(tree_file, "newick")
    except Exception as e:
        print(f"錯誤：讀取或解析樹檔案 '{tree_file}' 失敗: {e}")
        return

    print("開始為 PAML 準備輸入檔案...")
    
    gene_files = [f for f in os.listdir(aligned_cds_dir) if f.endswith(".fasta")]
    
    for gene_file in gene_files:
        gene_name = os.path.splitext(gene_file)[0]
        input_fasta_path = os.path.join(aligned_cds_dir, gene_file)
        gene_dir = os.path.join(paml_input_dir, gene_name)
        os.makedirs(gene_dir)
        output_phylip_path = os.path.join(gene_dir, f"{gene_name}.phy")
        try:
            # 讀取 FASTA 比對文件
            alignment = AlignIO.read(input_fasta_path, "fasta")
            
            # 【核心修改】移除了遍歷 alignment 並替換 record.id 的代碼塊
            # 現在直接使用原始的物種名

            # 使用優化後的函數寫入 PHYLIP 文件
            write_sequential_phylip(alignment, output_phylip_path)

        except Exception as e:
            print(f"  - 警告：處理基因 '{gene_name}' 時出錯: {e}，已跳過。")
            if os.path.exists(gene_dir): shutil.rmtree(gene_dir)
            continue
            
        # 將帶有原始物種名的樹文件寫入每個基因的目錄
        output_tree_path = os.path.join(gene_dir, "species_tree_original_names.nwk")
        write_paml_newick(tree, output_tree_path)
        
        # 生成 codeml 的控制文件 (.ctl)
        ctl_content = f"""
        seqfile = {gene_name}.phy
       treefile = species_tree_original_names.nwk
        outfile = {gene_name}_one_ratio.out

        noisy = 3
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 0
    fix_omega = 0
        omega = .4
        """
        ctl_path = os.path.join(gene_dir, "codeml_one_ratio.ctl")
        with open(ctl_path, "w") as f_ctl:
            # 清理字符串前後的空白並按行寫入
            f_ctl.write("".join(line.strip() + "\n" for line in ctl_content.strip().split("\n")))
            
        print(f"  -> 已成功為基因 '{gene_name}' 準備好所有輸入檔案。")

    # 【核心修改】移除了生成 "species_name_mapping.txt" 的代碼
    print(f"\nPAML 輸入檔案準備完成，所有文件均使用原始物種名。")

if __name__ == "__main__":
    prepare_paml_inputs_final_final()