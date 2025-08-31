import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Phylo
import os

def plot_final_heatmap():
    """
    主函數，讀取RSCU矩陣和系統發育樹，並繪製最終可調整細節的熱圖。
    """
    # --- 1. 定義輸入文件 ---
    rscu_csv_file = "rscu_matrix.csv"
    tree_file = "species_tree.nwk"
    output_image_file = "final_heatmap_customized.png"

    # --- 2. 【新增】定義圖形調整參數 ---
    # 您可以在這裡修改這些值來自訂您的圖片
    dendrogram_linewidth = 1.5  # 左側聚類樹線條的寬度，可以設為 1, 1.5, 2 等
    codon_label_fontsize = 15    # 右側密碼子標籤的字體大小，可以設為 6, 7, 8 等

    # --- 3. 檢查輸入文件是否存在 (與之前相同) ---
    if not os.path.exists(rscu_csv_file):
        print(f"錯誤: RSCU矩陣文件 '{rscu_csv_file}' 不存在。")
        return
    if not os.path.exists(tree_file):
        print(f"錯誤: 系統發育樹文件 '{tree_file}' 不存在。")
        return

    # --- 4. 讀取數據和系統發育樹 (與之前相同) ---
    print("正在讀取RSCU矩陣和系統發育樹...")
    rscu_df = pd.read_csv(rscu_csv_file, index_col=0)
    try:
        tree = Phylo.read(tree_file, "newick")
    except Exception as e:
        print(f"錯誤: 解析樹文件 '{tree_file}' 時出錯: {e}")
        return

    # --- 5. 根據樹的順序對RSCU矩陣的列進行重排 (與之前相同) ---
    tree_species_order = [leaf.name for leaf in tree.get_terminals()]
    matrix_species = list(rscu_df.columns)
    if sorted(tree_species_order) != sorted(matrix_species):
        print("錯誤: 樹文件中的物種名和RSCU矩陣中的物種名不完全匹配！")
        return
    
    rscu_df_ordered = rscu_df[tree_species_order]
    print("已根據系統發育樹的順序重排物種列。")

    # --- 6. 繪製熱圖 ---
    print("正在生成最終的熱圖...")
    
    g = sns.clustermap(
        rscu_df_ordered,
        col_cluster=False,
        row_cluster=True,
        cmap="viridis",
        linewidths=0.5,
        yticklabels=True,
        xticklabels=True,
        figsize=(18, 22),
        # **【核心修改1】** 透過 tree_kws 參數調整聚類樹線條寬度
        tree_kws={'linewidths': dendrogram_linewidth}
    )

    # **【核心修改2】** 調整Y軸（密碼子）標籤的字體大小
    # g.ax_heatmap 指的是中間熱圖所在的子圖 (Axes)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), size=codon_label_fontsize)
    
    # 旋轉X軸的標籤以防重疊
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    
    g.fig.suptitle("Codon Usage Heatmap Ordered by Phylogeny", fontsize=16)

    # --- 7. 保存最終的圖像 ---
    plt.savefig(output_image_file, dpi=300, bbox_inches='tight')
    print(f"\n分析完成！最終圖像已保存到 '{output_image_file}'。")


if __name__ == "__main__":
    plot_final_heatmap()