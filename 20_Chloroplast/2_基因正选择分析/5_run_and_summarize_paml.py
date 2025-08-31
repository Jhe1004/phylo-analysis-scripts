import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
import numpy as np

def generate_results_summary_english(results_df):
    """
    【英文版】根據分析結果生成一段可用於論文的、專業的英文描述性文字。
    """
    # 獲取基本信息
    num_genes = len(results_df['gene'].unique())
    num_lineages = len(results_df['foreground_label'].unique())
    
    # 篩選出 P<0.05 的顯著結果
    significant_df = results_df[results_df['p_value'] < 0.05].copy()
    
    # 從顯著結果中，再篩選出 Ka/Ks > 1 的正選擇結果
    positive_selection_df = significant_df[significant_df['omega_foreground'] > 1].copy()
    
    # --- 開始構建英文文本 ---
    summary_lines = []
    
    # --- 第一段：引言和總體趨勢 ---
    para1 = (f"To investigate the selective pressures on protein-coding genes within key clades, a series of branch-model "
             f"analyses were conducted using PAML on {num_genes} genes across {num_lineages} designated foreground lineages. "
             f"The results indicate that the majority of these genes are evolving under strong purifying selection (ω < 1) "
             f"across all tested lineages, consistent with their conserved functions.")
    summary_lines.append(para1)

    # --- 第二段：描述正選擇的結果 ---
    if not positive_selection_df.empty:
        para2_intro = (f"However, we identified {len(positive_selection_df)} instance(s) where genes exhibited "
                       f"significant evidence of positive selection (ω > 1, P < 0.05) in the foreground branches.")
        
        details = []
        # 使用斜體和更好的格式來描述基因
        positive_selection_df['gene_formatted'] = positive_selection_df['gene'].apply(lambda x: f"*{x}*")
        
        for index, row in positive_selection_df.iterrows():
            omega_str = ">> 2.0" if row['omega_foreground'] >= 2 else f"{row['omega_foreground']:.3f}"
            p_val_str = f"< 0.001" if row['p_value'] < 0.001 else f"{row['p_value']:.3f}"
            details.append(f"the gene {row['gene_formatted']} in the {row['foreground_label']} (ω = {omega_str}, P = {p_val_str})")
        
        # 根據細節數量選擇合適的連接詞
        if len(details) == 1:
            para2_details = f"Specifically, significant positive selection was detected for {details[0]}."
        else:
            details_string = "; ".join(details[:-1]) + f"; and {details[-1]}"
            para2_details = f"These instances include: {details_string}."
            
        summary_lines.append(para2_intro + " " + para2_details)
    else:
        para2 = "Across all tested foreground lineages, no genes were found to be under significant positive selection."
        summary_lines.append(para2)
        
    # --- 第三段：結論 ---
    para3 = ("These findings suggest that while most plastid genes are highly conserved, certain genes may have "
             "undergone episodic adaptive evolution in specific lineages.")
    summary_lines.append(para3)
    
    # 使用段落（雙換行符）連接所有部分
    return "\n\n".join(summary_lines)


def plot_heatmap_final():
    """
    主函數：【最终版】
    繪製無數字標註的熱圖，並在終端生成結果描述。
    """
    results_csv_file = "paml_final_results.csv"
    output_image_file = "kaks_heatmap_final.png"
    summary_text_file = "results_summary_english.txt" # 指定結果描述的輸出文件名

    print("="*50)
    print("PAML 分析腳本第二部分：生成最終結果圖和英文描述")
    print("="*50)

    # --- 1. 讀取並檢查數據 ---
    if not os.path.exists(results_csv_file):
        print(f"錯誤：找不到結果文件 '{results_csv_file}'。")
        return
    print(f"正在從 '{results_csv_file}' 讀取數據...")
    results_df = pd.read_csv(results_csv_file)
    # ... (数据检查部分代码省略) ...

    # --- 2. 數據整理 ---
    print("正在整理數據...")
    label_map = {}
    for set_id, group in results_df.groupby('foreground_set_id'):
        first_species_list = group['foreground_species'].iloc[0]
        first_species_name = first_species_list.split(',')[0]
        genus_name = re.split(r'[_ ]', first_species_name)[0]
        label_map[set_id] = f"{genus_name}_clade"
    results_df['foreground_label'] = results_df['foreground_set_id'].map(label_map)

    # --- 3. 生成結果描述文字並保存到文件 ---
    print(f"正在生成結果描述並保存到 '{summary_text_file}'...")
    results_summary = generate_results_summary_english(results_df)
    with open(summary_text_file, 'w', encoding='utf-8') as f:
        f.write(results_summary)
    print(" -> 描述文件已保存。")

    # --- 4. 準備繪圖數據 ---
    heatmap_data = results_df.pivot(index='foreground_label', columns='gene', values='omega_foreground')
    p_values = results_df.pivot(index='foreground_label', columns='gene', values='p_value')
    mask = p_values >= 0.05
    heatmap_data[heatmap_data > 2] = 2
    
    # --- 5. 繪製熱圖 ---
    print("正在繪製熱圖...")
    mean_omega_per_clade = heatmap_data.mean(axis=1).sort_values(ascending=False)
    heatmap_data_sorted = heatmap_data.loc[mean_omega_per_clade.index]
    mask_sorted = mask.loc[mean_omega_per_clade.index]
    
    fig_width = max(10, len(heatmap_data_sorted.columns) * 0.4)
    fig_height = max(6, len(heatmap_data_sorted.index) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # 第一層：畫背景
    sns.heatmap(heatmap_data_sorted, annot=False, cmap=['lightgrey'], cbar=False, ax=ax)
    
    # 第二層：疊加顯著結果
    ax = sns.heatmap(
        heatmap_data_sorted,
        mask=mask_sorted,
        annot=False,  # **核心修改**：設置為 False，不顯示數字
        cmap="coolwarm",
        linewidths=.5,
        vmin=0,
        vmax=2,
        ax=ax,
        cbar_kws={'label': 'Ka/Ks (ω) Ratio'}
    )

    ax.set_title("Significant Ka/Ks (ω) Ratios (P < 0.05)", fontsize=16, pad=20)
    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel("Foreground Lineage", fontsize=12)
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)

    plt.savefig(output_image_file, dpi=300, bbox_inches='tight')
    print(f"\n熱圖已儲存到 '{output_image_file}'。")
    print("\n可視化與結果描述生成全部完成！")

if __name__ == "__main__":
    plot_heatmap_final()