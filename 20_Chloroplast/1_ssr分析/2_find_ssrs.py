import os
import subprocess
import shutil

def run_misa_analysis():
    """
    主函數，自動化執行 MISA 工具來尋找 SSRs。
    【最終完美版】解決了 MISA 輸出路徑不正確的問題。
    """
    input_directory = "output_genome_fasta"
    output_directory = "output_ssr_results_misa"

    # MISA 腳本和設定檔的路徑
    misa_script_path = os.path.abspath("misa.py")
    misa_ini_path = os.path.abspath("misa.ini")

    # --- 1. 檢查環境 ---
    if not (os.path.exists(misa_script_path) and os.path.exists(misa_ini_path)):
        print(f"錯誤：找不到 '{misa_script_path}' 或 '{misa_ini_path}'。")
        return
    

    if not os.path.isdir(input_directory):
        print(f"錯誤：輸入文件夾 '{input_directory}' 不存在。")
        return

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    print("\n開始使用 MISA 腳本在所有基因組中搜索 SSRs...")


    # --- 3. 遍歷 Fasta 檔案並執行 MISA ---
    for filename in os.listdir(input_directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            species_name = os.path.splitext(filename)[0]
            source_fasta_path = os.path.join(input_directory, filename)
            
            print(f"  -> ------------------------------------")
            print(f"  -> 正在處理: {species_name}")

            # **【核心修改】**
            # 1. 將 Fasta 檔案複製到輸出目錄
            temp_fasta_path = os.path.join(output_directory, filename)
            shutil.copy(source_fasta_path, temp_fasta_path)
            
            # 2. 建立並執行 MISA 命令，這次是對複製過去的檔案進行操作
            #    因為 misa.pl 和 misa.ini 都在主目錄，我們直接提供它們的絕對路徑
            #    MISA 會在 temp_fasta_path 所在的目錄（也就是輸出目錄）生成結果
            command = ["python", misa_script_path, temp_fasta_path]
            
            try:
                # 這次我們不需要切換工作目錄(cwd)
                result = subprocess.run(
                    command, 
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                    check=True
                )
                print(f"    - MISA 分析成功。結果已生成在 '{output_directory}'。")

            except subprocess.CalledProcessError as e:
                print(f"    - 錯誤: MISA 執行失敗。")
                print(f"    - MISA 輸出: {e.stderr}")
            finally:
                # 3. 無論成功或失敗，都刪除臨時複製的 Fasta 檔案，保持輸出目錄乾淨
                if os.path.exists(temp_fasta_path):
                    os.remove(temp_fasta_path)

    print("\n所有物種的 SSR 搜索已完成！")


if __name__ == "__main__":
    run_misa_analysis()