import os
import shutil
import subprocess
import re
from scipy.stats import chi2
import pandas as pd
from Bio import Phylo
import multiprocessing # 引入多进程库

def parse_codeml_output_ultimate(outfile_path):
    """
    解析 codeml 的 .out 檔案。
    (此函數無需修改)
    """
    lnl = None
    omegas = {}
    try:
        with open(outfile_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        lnl_match = re.search(r"lnL\(.*?np:\s*\d+\):\s*([-\d\.]+)", content)
        if lnl_match:
            lnl = float(lnl_match.group(1))
        two_ratio_match = re.search(r"w \(dN/dS\) for branches:\s+([\d\.]+)\s+([\d\.]+)", content)
        if two_ratio_match:
            omegas[0] = float(two_ratio_match.group(1))
            omegas[1] = float(two_ratio_match.group(2))
        else:
            one_ratio_match = re.search(r"omega\s+\(dN/dS\)\s+=\s+([.\d]+)", content)
            if one_ratio_match:
                omegas[0] = float(one_ratio_match.group(1))
    except (FileNotFoundError, Exception) as e:
        # 在並行環境中，直接返回錯誤信息比打印更好
        return None, {}, f"解析文件 '{os.path.basename(outfile_path)}' 時發生錯誤: {e}"
    return lnl, omegas, None

def run_single_codeml_task(args):
    """
    【新增】執行單個 codeml 計算任務的工作函數，為多進程設計。
    """
    # 從傳入的元組中解包參數
    gene_name, set_id, fg_set, paml_input_dir, codeml_executable = args
    gene_dir_path = os.path.join(paml_input_dir, gene_name)
    
    try:
        tree_path = os.path.join(gene_dir_path, "species_tree_original_names.nwk")
        if not os.path.exists(tree_path):
            return "error", gene_name, set_id, f"找不到樹文件 '{os.path.basename(tree_path)}'"

        tree = Phylo.read(tree_path, "newick")
        tree_string = tree.format("newick").strip()
        
        for sp in fg_set:
            if f"{sp}:" in tree_string:
                tree_string = tree_string.replace(f"{sp}:", f"{sp} #1:")

        marked_tree_filename = f"species_tree_marked_{set_id}.nwk"
        marked_tree_path = os.path.join(gene_dir_path, marked_tree_filename)
        with open(marked_tree_path, 'w') as f:
            f.write(f" {len(tree.get_terminals())} 1\n")
            f.write(tree_string + "\n")

        output_filename = f"{gene_name}_two_ratio_{set_id}.out"
        ctl_filename = f"codeml_two_ratio_{set_id}.ctl"
        ctl_path = os.path.join(gene_dir_path, ctl_filename)

        ctl_content = (
            f"seqfile = {gene_name}.phy\n"
            f"treefile = {marked_tree_filename}\n"
            f"outfile = {output_filename}\n\n"
            "noisy = 0\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\n"
            "model = 2\nNSsites = 0\nfix_omega = 0\nomega = .4\n"
        )
        with open(ctl_path, "w") as f_ctl: f_ctl.write(ctl_content)

        command = [codeml_executable, ctl_filename]
        # 使用與 Python 3.6 兼容的 subprocess 參數
        subprocess.run(
            command,
            cwd=gene_dir_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=True
        )
        return "success", gene_name, set_id, "計算完成"
    except subprocess.CalledProcessError as e:
        error_message = f"codeml 執行失敗。\n錯誤信息: {e.stderr}"
        return "error", gene_name, set_id, error_message
    except Exception as e:
        return "error", gene_name, set_id, f"發生未知錯誤: {e}"

def run_paml_analysis_parallel():
    """
    主函數：【多進程修改版】
    """
    # --- 用戶配置區 ---
    # 【請您修改這裡】設置您希望使用的 CPU 核心數
    # os.cpu_count() 會自動獲取您機器的核心數，是一個不錯的默認值
    try:
        num_processes = os.cpu_count()
    except NotImplementedError:
        num_processes = 4 # 如果無法獲取CPU核心數，默認為4
        
    paml_input_dir = "paml_input"
    final_results_file = "paml_final_results.csv"
    foreground_list_file = "foreground_species_list.txt"

    print("="*50)
    print("PAML 分析腳本 (多進程版)：執行計算與生成結果表")
    print("="*50)

    # --- 步驟 1: 讀取前景列表和準備環境 (與之前相同) ---
    if not os.path.exists(foreground_list_file):
        print(f"錯誤：前景物種列表文件 '{foreground_list_file}' 不存在！")
        return
    
    all_foreground_sets = []
    with open(foreground_list_file, 'r') as f:
        for line in f:
            species = [s.strip() for s in line.strip().split(',') if s.strip()]
            if species: all_foreground_sets.append(species)
    
    if not all_foreground_sets:
        print(f"錯誤：未能從 '{foreground_list_file}' 中讀取到任何前景物種組合。")
        return

    print(f"成功從 '{foreground_list_file}' 讀取到 {len(all_foreground_sets)} 组前景分支待測試。")

    codeml_executable = shutil.which("codeml")
    if not codeml_executable:
        print("錯誤：在您的系統環境變數中找不到 'codeml' 可執行文件。")
        return

    if not os.path.isdir(paml_input_dir):
        print(f"錯誤：找不到 PAML 輸入目錄 '{paml_input_dir}'。")
        return
        
    gene_dirs = sorted([d for d in os.listdir(paml_input_dir) if os.path.isdir(os.path.join(paml_input_dir, d))])
    if not gene_dirs:
        print(f"警告：在 '{paml_input_dir}' 目錄下沒有找到任何基因文件夾。")
        return

    # --- 步驟 2: 創建所有待執行的任務列表 ---
    tasks = []
    for gene_name in gene_dirs:
        for i, fg_set in enumerate(all_foreground_sets):
            set_id = f"set_{i}"
            output_filepath = os.path.join(paml_input_dir, gene_name, f"{gene_name}_two_ratio_{set_id}.out")
            
            # 如果結果已存在，則不加入任務列表
            if not os.path.exists(output_filepath):
                task_args = (gene_name, set_id, fg_set, paml_input_dir, codeml_executable)
                tasks.append(task_args)

    if not tasks:
        print("\n所有計算任務的結果均已存在，無需執行新的計算。")
    else:
        print(f"\n檢測到 {len(tasks)} 個新的計算任務。")
        print(f"將使用 {num_processes} 個進程並行執行...")
        
        # --- 步驟 3: 使用多進程池執行所有任務 ---
        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.map(run_single_codeml_task, tasks)

        # --- 步驟 4: 報告執行結果 ---
        print("\n並行計算完成。結果摘要：")
        success_count = 0
        error_count = 0
        for status, gene_name, set_id, message in results:
            if status == "success":
                success_count += 1
            else:
                error_count += 1
                print(f"  - 失敗: 基因 '{gene_name}' - {set_id}。原因: {message}")
        print(f"  -> 總計：{success_count} 個任務成功，{error_count} 個任務失敗。")


    # --- 步驟 5: 順序解析所有結果並進行 LRT 檢定 ---
    print("\n開始解析所有結果並進行 LRT 檢定...")
    results_list = []
    for gene_name in gene_dirs:
        for i, fg_set in enumerate(all_foreground_sets):
            set_id = f"set_{i}"
            one_ratio_out = os.path.join(paml_input_dir, gene_name, f"{gene_name}_one_ratio.out")
            two_ratio_out = os.path.join(paml_input_dir, gene_name, f"{gene_name}_two_ratio_{set_id}.out")
            
            if not os.path.exists(one_ratio_out) or not os.path.exists(two_ratio_out):
                continue

            lnL0, omegas0, err0 = parse_codeml_output_ultimate(one_ratio_out)
            lnL1, omegas1, err1 = parse_codeml_output_ultimate(two_ratio_out)

            if lnL0 is None or lnL1 is None or not omegas0 or len(omegas1) < 2:
                continue
                
            lrt_stat = 2 * (lnL1 - lnL0)
            p_value = chi2.sf(lrt_stat, df=1) if lrt_stat > 0 else 1.0
            
            results_list.append({
                "gene": gene_name,
                "foreground_set_id": set_id,
                "foreground_species": ",".join(fg_set),
                "lnL_null": lnL0,
                "lnL_alt": lnL1,
                "omega_null": omegas0.get(0),
                "omega_foreground": omegas1.get(1),
                "omega_background": omegas1.get(0),
                "LRT_statistic": lrt_stat,
                "p_value": p_value
            })
    
    if not results_list:
        print("\n錯誤：未能成功解析任何基因的結果，無法生成報告。")
        return

    results_df = pd.DataFrame(results_list)
    results_df.to_csv(final_results_file, index=False, float_format="%.6f")
    print(f"\n分析完成！最終統計結果已儲存到 '{final_results_file}'。")
    print("您現在可以運行第二個腳本來繪製熱圖。")

if __name__ == "__main__":
    run_paml_analysis_parallel()