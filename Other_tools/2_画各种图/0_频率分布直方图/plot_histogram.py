#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/2_画各种图/0_频率分布直方图/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/2_画各种图/0_频率分布直方图/output"
INPUT_CSV_FILE = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/Other_tools/2_画各种图/0_频率分布直方图/input/data.csv"
CHART_TITLE = "频率分布直方图"
X_LABEL = "数值"
Y_LABEL = "频数"


OUTPUT_IMAGE_FILE = "distribution_plot.png"


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    input_file = INPUT_DIR / INPUT_CSV_FILE
    if not input_file.exists():
        raise FileNotFoundError(f"未找到输入 CSV 文件: {input_file}")
    df = pd.read_csv(input_file, header=None)
    data = df.iloc[:, -1].dropna()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(12, 7))
    sns.set_theme(style="whitegrid", palette="viridis")
    sns.histplot(data, bins="auto", kde=True, color="#5B4E82", line_kws={"linewidth": 2.5})
    plt.title(CHART_TITLE, fontsize=16)
    plt.xlabel(X_LABEL, fontsize=12)
    plt.ylabel(Y_LABEL, fontsize=12)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / OUTPUT_IMAGE_FILE, dpi=300)
    plt.close()
    print("直方图绘制完成。")


if __name__ == "__main__":
    main()
