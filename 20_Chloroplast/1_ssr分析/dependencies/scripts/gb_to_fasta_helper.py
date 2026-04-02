import os
import re
from Bio import SeqIO

def clean_name(name):
    """
    清理名稱字符串，將空格和特殊字符替換為下划線。
    """
    special_chars = r" |\(|\)|\[|\]"
    cleaned = re.sub(special_chars, "_", name)
    cleaned = re.sub(r"__+", "_", cleaned)
    cleaned = cleaned.strip('_')
    return cleaned

def convert_gb_to_fasta_robust():
    """
    將GenBank文件轉換為Fasta格式的全基因組序列文件。
    【格式強化版】確保輸出的 Fasta 格式最標準、最純淨，以避免 MISA 解析錯誤。
    """
    input_directory = "input_gb_files"
    output_directory = "output_genome_fasta"

    if not os.path.isdir(input_directory):
        print(f"錯誤：輸入文件夾 '{input_directory}' 不存在。")
        return

    if os.path.exists(output_directory):
        # 如果輸出文件夾已存在，先清空，確保是全新的開始
        for f in os.listdir(output_directory):
            os.remove(os.path.join(output_directory, f))
        print(f"已清空舊的輸出文件夾 '{output_directory}'。")
    else:
        os.makedirs(output_directory)
        print(f"已創建輸出文件夾：'{output_directory}'")

    print("\n開始將 GenBank 轉換為標準化 Fasta 格式...")
    converted_files = 0
    for filename in os.listdir(input_directory):
        if filename.endswith(".gb") or filename.endswith(".gbk"):
            input_filepath = os.path.join(input_directory, filename)
            
            original_base_name = os.path.splitext(filename)[0]
            species_name = clean_name(original_base_name)
            
            output_filename = f"{species_name}.fasta"
            output_filepath = os.path.join(output_directory, output_filename)

            try:
                # 讀取 GenBank 檔案
                record = SeqIO.read(input_filepath, "genbank")
                
                # **【核心修改】**
                # 1. 獲取純淨的序列字串
                sequence_str = str(record.seq).upper()
                # 2. 過濾掉所有非 ATCGN 的字元
                cleaned_sequence = re.sub(r'[^ATCGN]', '', sequence_str)
                
                # 3. 手動寫入標準 Fasta 格式
                with open(output_filepath, "w") as f_out:
                    # 寫入非常簡潔的標頭
                    f_out.write(f">{species_name}\n")
                    # 將序列每 60 個字元換行一次
                    for i in range(0, len(cleaned_sequence), 60):
                        f_out.write(cleaned_sequence[i:i+60] + "\n")
                
                print(f"  -> 成功轉換並標準化: '{filename}' -> '{output_filename}'")
                converted_files += 1

            except Exception as e:
                print(f"  - 錯誤: 處理檔案 {filename} 時發生錯誤: {e}")

    print(f"\n轉換完成！總共成功生成了 {converted_files} 個標準化的 Fasta 文件。")


if __name__ == "__main__":
    convert_gb_to_fasta_robust()