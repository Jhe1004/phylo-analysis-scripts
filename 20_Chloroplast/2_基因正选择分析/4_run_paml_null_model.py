import os
import subprocess
import shutil

def run_codeml_null_model():
    """
    主函數，自動化地為所有基因執行 codeml 的單比例模型（虛無假設）。
    """
    paml_input_dir = "paml_input"
    
    # --- 1. 檢查環境和目錄 ---
    # 使用 shutil.which 來自動在環境路徑中尋找 codeml
    codeml_executable = shutil.which("codeml")
    if not codeml_executable:
        print("錯誤：在您的系統環境中找不到 'codeml' 可執行程式。")
        print("請確認您已在 conda 環境中成功安裝 PAML (conda install -c bioconda paml)。")
        return
        
    if not os.path.isdir(paml_input_dir):
        print(f"錯誤：找不到 PAML 輸入文件夾 '{paml_input_dir}'。請先執行 B2_prep_paml_files.py。")
        return

    print("開始執行 codeml 進行選擇壓力分析 (虛無模型：one-ratio)...")
    print(f"找到 codeml 執行檔位於: {codeml_executable}")
    print("這個過程可能會非常耗時，請耐心等待...")

    # --- 2. 遍歷所有基因的獨立文件夾 ---
    gene_dirs = [d for d in os.listdir(paml_input_dir) if os.path.isdir(os.path.join(paml_input_dir, d))]
    total_genes = len(gene_dirs)
    
    for i, gene_name in enumerate(sorted(gene_dirs), 1):
        gene_dir_path = os.path.join(paml_input_dir, gene_name)
        control_file = "codeml_one_ratio.ctl"
        control_file_path = os.path.join(gene_dir_path, control_file)
        
        if not os.path.exists(control_file_path):
            print(f"  - 警告：在 '{gene_dir_path}' 中找不到控制檔 '{control_file}'，已跳過。")
            continue
            
        print(f"\n--- 正在處理基因: {gene_name} ({i}/{total_genes}) ---")
        
        # 建立 codeml 命令
        # codeml 會讀取它所在目錄下的控制檔，所以我們只需要傳遞控制檔名
        command = [codeml_executable, control_file]
        
        try:
            # 關鍵：我們必須在該基因自己的文件夾中執行 codeml 命令
            # 這樣 codeml 才能根據控制檔中的相對路徑找到 .phy 和 .nwk 檔案
            result = subprocess.run(
                command,
                cwd=gene_dir_path, # 設定命令的執行目錄
                capture_output=True, # 捕獲輸出
                text=True,           # 將輸出解碼為文字
                check=True           # 如果命令出錯則拋出例外
            )
            # codeml 的輸出很多，通常我們只關心它是否成功
            # 如果您想看詳細日誌，可以取消下面這行的註解
            # print(result.stdout)
            print(f"  -> 成功完成基因 '{gene_name}' 的計算。")

        except subprocess.CalledProcessError as e:
            print(f"  - 錯誤：在處理基因 '{gene_name}' 時 codeml 執行失敗。")
            print(f"  - PAML 錯誤訊息:\n{e.stderr}")
        # 處理與 Python 3.6 的相容性
        except TypeError:
            try:
                result = subprocess.run(
                    command,
                    cwd=gene_dir_path,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                    check=True
                )
                print(f"  -> 成功完成基因 '{gene_name}' 的計算 (Python 3.6 相容模式)。")
            except subprocess.CalledProcessError as e_compat:
                print(f"  - 錯誤：在處理基因 '{gene_name}' 時 codeml 執行失敗 (Python 3.6 相容模式)。")
                print(f"  - PAML 錯誤訊息:\n{e_compat.stderr}")

    print("\n所有基因的虛無模型計算已完成！")

if __name__ == "__main__":
    run_codeml_null_model()