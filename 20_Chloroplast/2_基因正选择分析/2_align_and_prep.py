import os
import subprocess
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def back_translate_alignment(aligned_pep_file, unaligned_cds_file):
    """
    使用比對好的蛋白質序列作為模板，來生成密碼子比對。
    這是 PAL2NAL 工具的核心邏輯。
    """
    aligned_pep = {rec.id: str(rec.seq) for rec in SeqIO.parse(aligned_pep_file, "fasta")}
    unaligned_cds = {rec.id: str(rec.seq) for rec in SeqIO.parse(unaligned_cds_file, "fasta")}
    
    aligned_cds_records = []
    
    for species_id, pep_seq in aligned_pep.items():
        if species_id not in unaligned_cds:
            continue
            
        cds_seq = unaligned_cds[species_id]
        aligned_cds_seq = ""
        cds_pos = 0
        
        for amino_acid in pep_seq:
            if amino_acid == '-':
                aligned_cds_seq += "---"
            else:
                codon = cds_seq[cds_pos : cds_pos + 3]
                aligned_cds_seq += codon
                cds_pos += 3
        
        new_record = SeqRecord(Seq(aligned_cds_seq), id=species_id, description="")
        aligned_cds_records.append(new_record)
        
    return aligned_cds_records

def align_and_prepare():
    """
    主函數，使用 MAFFT 進行蛋白質序列比對，然後用 Python 進行密碼子比對，
    並生成 ParaAT 所需的輸入檔案。
    【Python 3.6 相容最終版】
    """
    cds_dir = "output_cds_genes"
    pep_dir = "output_pep_genes"
    aligned_cds_dir = "aligned_cds"
    aligned_pep_dir = "aligned_pep"
    para_at_input_file = "homologs.txt"

    if shutil.which("mafft") is None:
        print("錯誤：找不到 'mafft' 命令。")
        return

    for d in [cds_dir, pep_dir]:
        if not os.path.isdir(d):
            print(f"錯誤：找不到輸入文件夾 '{d}'。")
            return
            
    for d in [aligned_cds_dir, aligned_pep_dir]:
        if not os.path.exists(d):
            os.makedirs(d)

    print("開始進行序列比對...")
    
    gene_files = [f for f in os.listdir(cds_dir) if f.endswith(".fasta")]
    
    for gene_file in gene_files:
        gene_name = os.path.splitext(gene_file)[0]
        print(f"  -> 正在處理基因: {gene_name}")

        input_pep_path = os.path.join(pep_dir, gene_file)
        aligned_pep_path = os.path.join(aligned_pep_dir, gene_file)
        input_cds_path = os.path.join(cds_dir, gene_file)
        aligned_cds_path = os.path.join(aligned_cds_dir, gene_file)

        try:
            pep_align_command = ["mafft", "--auto", input_pep_path]
            with open(aligned_pep_path, "w") as f_out:
                # **【核心修改】**
                # 將 subprocess.run 的參數替換為 Python 3.6 相容的寫法
                subprocess.run(
                    pep_align_command, 
                    stdout=f_out, 
                    stderr=subprocess.PIPE, 
                    check=True, 
                    universal_newlines=True # 使用 universal_newlines 替代 text=True
                )

            aligned_cds_records = back_translate_alignment(aligned_pep_path, input_cds_path)
            SeqIO.write(aligned_cds_records, aligned_cds_path, "fasta")

        except subprocess.CalledProcessError as e:
            print(f"    - 錯誤: 在比對基因 '{gene_name}' 的蛋白質時 MAFFT 執行失敗。")
            print(f"    - 錯誤訊息: {e.stderr}")
            continue
        except Exception as e:
            print(f"    - 發生未知錯誤: {e}")
            continue

    print("\n比對完成。正在生成 ParaAT 的輸入檔案 'homologs.txt'...")
    
    aligned_gene_files = [f for f in os.listdir(aligned_cds_dir) if f.endswith(".fasta")]

    try:
        with open(para_at_input_file, "w") as f_out:
            f_out.write("gene_name\tcds_path\tpep_path\n")
            for gene_file in aligned_gene_files:
                gene_name = os.path.splitext(gene_file)[0]
                cds_path = os.path.join(aligned_cds_dir, gene_file)
                pep_path = os.path.join(aligned_pep_dir, gene_file)
                f_out.write(f"{gene_name}\t{cds_path}\t{pep_path}\n")
                
        print(f"成功生成 '{para_at_input_file}'。")
        print("\n第二步已完成！")

    except Exception as e:
        print(f"錯誤：生成 '{para_at_input_file}' 時發生錯誤: {e}")

if __name__ == "__main__":
    align_and_prepare()