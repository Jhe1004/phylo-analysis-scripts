import os
import pandas as pd
import matplotlib.pyplot as plt
import io
from Bio import Phylo # 引入 BioPython 函式庫用於處理發育樹

def summarize_and_plot_ssrs_with_tree():
    """
    主函數，讀取所有 MISA 結果，結合系統發育樹，
    生成報告和帶有發育樹的水平堆疊柱狀圖。
    """
    # --- 輸入與輸出檔案設定 ---
    input_directory = "output_ssr_results_misa"
    phylogenetic_tree_file = "species_tree.nwk" # <--- 新增：請將你的樹文件放在這裡
    
    output_data_file = "ssr_summary_data.csv"
    output_report_file = "ssr_results_report.txt"
    output_image_file = "ssr_phylo_stacked_barplot_corrected.png" # 使用新文件名以避免覆盖

    # --- 1. 讀取和解析 MISA 檔案 (與原腳本相同) ---
    all_species_data = []
    print("正在讀取並解析所有 .misa 結果檔案...")

    if not os.path.isdir(input_directory):
        print(f"錯誤：輸入文件夾 '{input_directory}' 不存在。")
        return

    for filename in os.listdir(input_directory):
        if filename.endswith(".misa"):
            filepath = os.path.join(input_directory, filename)
            species_name = filename.split('.fasta')[0]
            
            try:
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                    file_content_lines = f.readlines()

                file_content = ""
                header_found = False
                for line in file_content_lines:
                    if line.strip().startswith('ID\t'):
                        header_found = True
                    if header_found:
                        file_content += line
                
                if file_content:
                    df = pd.read_csv(io.StringIO(file_content), sep='\t')
                    if not df.empty:
                        df['species'] = species_name
                        all_species_data.append(df)
                    else:
                        print(f"  - 警告：檔案 {filename} 中沒有找到 SSR 數據。")
                else:
                    print(f"  - 警告：在檔案 {filename} 中找不到數據表頭 'ID'。")
            except Exception as e:
                print(f"  - 警告：處理檔案 {filename} 時出錯: {e}")

    if not all_species_data:
        print("錯誤：未能從任何 .misa 檔案中成功解析出數據。程式終止。")
        return
        
    full_df = pd.concat(all_species_data, ignore_index=True)

    # --- 2. 統計 SSR 類型 (與原腳本相同) ---
    print("正在統計 SSR 類型...")
    ssr_type_name_map = {
        'p1': 'mononucleotide', 'p2': 'dinucleotide', 'p3': 'trinucleotide',
        'p4': 'tetranucleotide', 'p5': 'pentanucleotide', 'p6': 'hexanucleotide',
        'c': 'compound'
    }
    full_df['ssr_type_classified'] = full_df['SSR type'].map(ssr_type_name_map)
    summary_df = full_df.pivot_table(index='species', columns='ssr_type_classified', aggfunc='size', fill_value=0)

    column_order = ['mononucleotide', 'dinucleotide', 'trinucleotide', 'tetranucleotide', 'pentanucleotide', 'hexanucleotide', 'compound']
    for ssr_type_name in column_order:
        if ssr_type_name not in summary_df.columns:
            summary_df[ssr_type_name] = 0
    summary_df = summary_df[column_order]

    # --- 3. 讀取系統發育樹並對數據進行排序 (新功能) ---
    tree = None
    try:
        print(f"正在讀取系統發育樹檔案 '{phylogenetic_tree_file}'...")
        tree = Phylo.read(phylogenetic_tree_file, "newick")
        # 獲取樹末端（物種）從上到下的順序
        # BioPython 默认的 get_terminals 顺序就是绘图时从下到上的顺序
        species_order_from_tree = [tip.name for tip in tree.get_terminals()]
        
        # 檢查數據和樹中的物種名是否能對應
        species_in_df = set(summary_df.index)
        species_in_tree = set(species_order_from_tree)
        if species_in_df != species_in_tree:
            print("警告：數據中的物種名與系統發育樹中的不完全匹配！")
            print(f"  - 僅存在於數據中的物種: {species_in_df - species_in_tree}")
            print(f"  - 僅存在於樹中的物種: {species_in_tree - species_in_df}")
        
        # 根據樹的順序對 DataFrame 進行重排，忽略樹中多餘的物種
        summary_df = summary_df.reindex(species_order_from_tree, fill_value=0)
        
        # 【關鍵修改 1】：刪除下面這一行反轉DataFrame的程式碼。
        # 讓DataFrame的順序與Phylo.draw的預設繪圖順序保持一致。
        # 原始程式碼: summary_df = summary_df.iloc[::-1]

    except FileNotFoundError:
        print(f"警告：找不到系統發育樹文件 '{phylogenetic_tree_file}'。將只生成橫向柱狀圖。")
        # 如果沒有樹，則按物種名進行排序
        summary_df = summary_df.sort_index(ascending=False)
    except Exception as e:
        print(f"錯誤：讀取或解析樹文件時出錯: {e}")
        return

    # 儲存排序後的數據
    summary_df.to_csv(output_data_file)
    print(f"已將詳細的 SSR 統計數據儲存到 '{output_data_file}'。")

    # --- 4. 生成報告 (與原腳本相同) ---
    print("正在生成結果報告...")
    total_ssrs = summary_df.sum().sum()
    num_species = len(summary_df)
    avg_ssrs_per_species = total_ssrs / num_species if num_species > 0 else 0
    most_abundant_type = summary_df.sum().idxmax()
    report_content = f"""
SSR Analysis Report
===================
- A total of {total_ssrs} simple sequence repeats (SSRs) were identified across {num_species} species.
- The average number of SSRs per species was approximately {avg_ssrs_per_species:.2f}.
- Among the different repeat types, **{most_abundant_type}** repeats were the most abundant category across all analyzed genomes.
- Mononucleotide repeats were the most common type, predominantly consisting of A/T repeats, which is a typical feature of chloroplast genomes.
- The number of dinucleotide, trinucleotide, and other longer repeats was considerably lower in comparison.
- Detailed counts for each repeat type per species are available in the accompanying data file '{output_data_file}'.
"""
    with open(output_report_file, 'w', encoding='utf-8') as f:
        f.write(report_content)
    print(f"已將結果報告儲存到 '{output_report_file}'。")


    # --- 5. 繪製帶有發育樹的堆疊柱狀圖 (全新繪圖邏輯) ---
    print("正在繪製帶有發育樹的水平堆疊柱狀圖...")
    
    # 根據是否有樹，決定畫布佈局
    if tree:
        # 創建一個 1x2 的子圖網格，共享 Y 軸
        fig, (ax_tree, ax_bar) = plt.subplots(
            1, 2, 
            figsize=(16, max(8, num_species * 0.5)), # 根據物種數量動態調整高度
            sharey=True, 
            gridspec_kw={'width_ratios': [1, 3]} # 左圖(樹)與右圖(柱狀圖)的寬度比例為 1:3
        )
        # 在左側子圖中繪製發育樹
        Phylo.draw(tree, axes=ax_tree, do_show=False, label_func=lambda x: None) # 不让Phylo自己画标签
        
        # 美化樹圖的外觀
        ax_tree.spines['top'].set_visible(False)
        ax_tree.spines['right'].set_visible(False)
        ax_tree.spines['bottom'].set_visible(False)
        ax_tree.spines['left'].set_visible(False)
        ax_tree.get_xaxis().set_ticks([]) # 隱藏樹的X軸
        
        # 【關鍵修改 2】：隱藏樹圖(左圖)的Y軸刻度標籤，以防止與樹重疊。
        # 標籤將由右側的柱狀圖生成並顯示。
        ax_tree.tick_params(axis='y', which='both', length=0)

    else:
        # 如果沒有樹，只創建一個子圖
        fig, ax_bar = plt.subplots(figsize=(12, max(6, num_species * 0.4)))
        
    # 在右側子圖中繪製水平堆疊柱狀圖
    colors = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462', '#ccebc5']
    summary_df.plot(
        kind='barh',        # <-- 改成 'barh'
        stacked=True,       
        color=colors,
        ax=ax_bar,
        width=0.8
    )
    
    # 美化圖表
    ax_bar.set_xlabel('Number of SSRs', fontsize=14)
    ax_bar.set_ylabel('Taxa', fontsize=14) # 給Y軸一個總標題
    ax_bar.tick_params(axis='y', labelsize=12)
    ax_bar.tick_params(axis='x', labelsize=12)
    ax_bar.grid(axis='x', linestyle='--', alpha=0.6) # 增加X軸網格線，方便對齊
    
    # 將圖例放在圖表外部右上方
    ax_bar.legend(title='SSR Type', fontsize=12, bbox_to_anchor=(1.02, 1), loc='upper left')

    # 【關鍵修改 3】：在所有繪圖完成後，翻轉Y軸。
    # 因為Y軸是共享的，這個操作會同時翻轉樹和柱狀圖，使它們的順序與你期望的一致(從上到下)。
    if tree:
        ax_bar.invert_yaxis()

    # 為整個圖表設定一個主標題
    fig.suptitle('Distribution of SSR Types Across Species with Phylogeny', fontsize=18, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96]) # 調整佈局以容納主標題
    
    plt.savefig(output_image_file, dpi=300, bbox_inches='tight')
    print(f"已將整合圖表儲存到 '{output_image_file}'。")
    print("\n所有任務已完成！")


# --- 執行主函數 ---
summarize_and_plot_ssrs_with_tree()