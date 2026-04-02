#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
INPUT_CSV_FILE = "input.csv"
OUTPUT_IMAGE_FILE = "violin_plot_with_points.png"
PLOT_TITLE = "数据分布小提琴图"
Y_AXIS_LABEL = "数值"
PLOT_COLOR = "cornflowerblue"
SHOW_INDIVIDUAL_POINTS = True
POINT_COLOR = "black"
POINT_SIZE = 4
POINT_ALPHA = 0.6
SET_Y_AXIS_RANGE = True
Y_AXIS_MIN = 0.0
Y_AXIS_MAX = 1.0


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / INPUT_DIRECTORY
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def main() -> None:
    input_file = INPUT_DIR / INPUT_CSV_FILE
    if not input_file.exists():
        raise FileNotFoundError(f"未找到输入 CSV 文件: {input_file}")
    data = pd.read_csv(input_file)
    plot_data = data.iloc[:, 0]
    column_name = data.columns[0]
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(8, 10))
    ax = sns.violinplot(y=plot_data, color=PLOT_COLOR)
    if SHOW_INDIVIDUAL_POINTS:
        sns.stripplot(y=plot_data, color=POINT_COLOR, s=POINT_SIZE, alpha=POINT_ALPHA, jitter=0.1, ax=ax)
    ax.set_title(PLOT_TITLE, fontsize=16)
    ax.set_ylabel(Y_AXIS_LABEL, fontsize=12)
    ax.set_xlabel(column_name, fontsize=12)
    if SET_Y_AXIS_RANGE:
        ax.set_ylim(Y_AXIS_MIN, Y_AXIS_MAX)
    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT_DIR / OUTPUT_IMAGE_FILE, dpi=300)
    plt.close()
    print("小提琴图绘制完成。")


if __name__ == "__main__":
    main()
