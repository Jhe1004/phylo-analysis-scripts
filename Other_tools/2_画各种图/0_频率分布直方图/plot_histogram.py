# plot_histogram_simple_fixed.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import urllib.request
from matplotlib.font_manager import FontProperties

# --- 用户配置区域 ---
# 请在这里修改您要处理的文件名和希望输出的图片名

# 1. 输入的CSV文件名
# 例如: "Node_1_gamma_diff.csv" 或 "gamma_differences_log.csv"
INPUT_CSV_FILE = "aa.csv"

# 2. 输出的图片文件名
# 您可以根据输入文件名来命名，例如: "distribution_plot_for_Node_1.png"
OUTPUT_IMAGE_FILE = "distribution_plot.png"

# 3. 图表标题 (可选, 您可以修改成任何您想要的文字)
CHART_TITLE = f"Gamma值差异频率分布\n(数据来源: {os.path.basename(INPUT_CSV_FILE)})"

# --- 配置结束 ---


def get_cjk_font():
    """
    下载并加载一个支持中日韩文字的开源字体。
    这里使用 Google Noto Sans SC (简体中文)。
    """
    font_url = 'https://raw.githubusercontent.com/google/fonts/main/ofl/notosanssc/NotoSansSC-Regular.otf'
    font_filename = "NotoSansSC-Regular.otf"

    if not os.path.exists(font_filename):
        try:
            print(f"正在下载中文字体: {font_filename}...")
            urllib.request.urlretrieve(font_url, font_filename)
            print("字体下载完成。")
        except Exception as e:
            print(f"错误：无法下载字体。请检查网络连接。错误信息: {e}")
            # 如果下载失败，返回None，后续将使用默认字体
            return None
            
    # 加载下载好的字体文件
    return FontProperties(fname=font_filename)


def create_histogram(
    input_file: str,
    output_file: str,
    title: str,
    bins: int = None,
    ylabel: str = "频数 (Frequency)",
    xlabel: str = '数值 (Value)' # 为X轴也定义一个中文标签
):
    """
    读取CSV文件并自动为最后一列数据绘制频率分布直方图。

    参数:
    - input_file (str): 输入的CSV文件名。
    - output_file (str): 输出的图片文件名。
    - title (str): 图表标题。
    - bins (int, optional): 直方图的分箱数。如果为None，则自动确定。
    - ylabel (str, optional): Y轴标签。
    - xlabel (str, optional): X轴标签。
    """
    # --- 0. 加载中文字体 ---
    cjk_font = get_cjk_font()
    if cjk_font is None:
        print("警告: 未能加载中文字体，图表中的中文可能无法正常显示。")

    # --- 1. 数据加载与验证 ---
    try:
        df = pd.read_csv(input_file, header=None)
        column_name = df.columns[-1]
        print(f"自动选择最后一列 (列索引: {column_name}) 进行分析")

        data = df[column_name].dropna()
        if data.empty:
            print(f"错误: 列 '{column_name}' 中没有有效数据可供绘制。")
            return
    except FileNotFoundError:
        print(f"错误: 文件 '{input_file}' 未找到。请检查上面的 INPUT_CSV_FILE 变量是否正确。")
        return
    except IndexError:
        print(f"错误: CSV文件 '{input_file}' 为空或无法读取列。")
        return
    except Exception as e:
        print(f"读取或处理文件时出错: {e}")
        return

    # --- 2. 绘图设置 ---
    sns.set_theme(style="whitegrid", palette="viridis") 
    plt.figure(figsize=(12, 7)) 

    # --- 3. 绘制直方图和密度曲线 ---
    plot_bins = "auto" if bins is None else bins
    # 使用较深的颜色以提高可见性
    sns.histplot(data, bins=plot_bins, kde=True, color="#5B4E82", line_kws={'linewidth': 2.5})

    # --- 4. 设置标题和标签 ---
    # 【重要改动】在这里应用下载好的字体
    plt.title(title, fontsize=18, pad=20, fontproperties=cjk_font)
    plt.xlabel(xlabel, fontsize=14, fontproperties=cjk_font)
    plt.ylabel(ylabel, fontsize=14, fontproperties=cjk_font)
    
    # 解决使用自定义字体后，坐标轴负号可能显示异常的问题
    plt.rcParams['axes.unicode_minus'] = False 

    # 调整刻度标签的大小
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # 确保布局紧凑
    plt.tight_layout()

    # --- 5. 保存图片 ---
    try:
        plt.savefig(output_file, dpi=300)
        print(f"直方图已成功保存到: {output_file}")
    except Exception as e:
        print(f"保存图片时出错: {e}")

    plt.close()

if __name__ == "__main__":
    # 直接调用主函数，使用在脚本顶部定义的配置
    # 【重要改动】为X轴也传入一个中文标签，保持统一
    final_xlabel = "Gamma 值" 
    create_histogram(
        input_file=INPUT_CSV_FILE,
        output_file=OUTPUT_IMAGE_FILE,
        title=CHART_TITLE,
        xlabel=final_xlabel
    )