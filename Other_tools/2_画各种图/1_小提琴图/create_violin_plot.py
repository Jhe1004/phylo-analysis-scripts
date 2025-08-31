import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================================================================
#  参数设置 (User-configurable parameters)
# ==============================================================================

# 1. 输入和输出文件设置
# ------------------------------------------------------------------------------
input_file = 'input.csv'
output_file = 'violin_plot_with_points.png'


# 2. 图表样式设置
# ------------------------------------------------------------------------------
plot_title = '数据分布小提琴图 (含数据点)'
y_axis_label = '数值'
plot_color = 'cornflowerblue'
figure_width = 8
figure_height = 10 # 适当增加了高度，以便更好地展示数据点和标签


# 3. 数据点显示设置 (新增功能)
# ------------------------------------------------------------------------------
# show_individual_points: 是否在图上显示所有的数据点。
# 设置为 True 来显示，设置为 False 来隐藏。
show_individual_points = True

# point_color: 数据点的颜色。
point_color = 'black'

# point_size: 数据点的大小。
point_size = 4

# point_alpha: 数据点的透明度（0.0-1.0）。
# 当数据点密集时，使用较低的透明度（如0.5）可以防止点重叠成一团黑。
point_alpha = 0.6


# 4. 坐标轴范围设置 (新增功能)
# ------------------------------------------------------------------------------
# set_y_axis_range: 是否自定义Y轴的范围。
# 设置为 True 来使用下面的min/max值，设置为 False 则由程序自动确定范围。
set_y_axis_range = True

# y_axis_min: Y轴的最小值。仅在 set_y_axis_range 为 True 时生效。
y_axis_min = 0.0

# y_axis_max: Y轴的最大值。仅在 set_y_axis_range 为 True 时生效。
y_axis_max = 1.0


# ==============================================================================
#  主程序 (Main script logic)
#  通常您不需要修改以下部分
# ==============================================================================

def create_violin_plot_v2():
    """
    读取CSV文件并生成带有数据点和可自定义Y轴范围的小提琴图。
    """
    try:
        # 读取数据
        print(f"正在读取数据文件: {input_file}...")
        data = pd.read_csv(input_file)
        
        if data.empty:
            print("错误：数据文件为空，无法生成图表。")
            return
            
        plot_data = data.iloc[:, 0]
        column_name = data.columns[0]
        print("数据读取成功，开始生成图表...")

        # 设置图表风格
        sns.set_theme(style="whitegrid")

        # 创建图表
        plt.figure(figsize=(figure_width, figure_height))
        
        # 1. 绘制小提琴图
        ax = sns.violinplot(y=plot_data, color=plot_color)

        # 2. (新增) 在小提琴图上叠加显示所有数据点
        if show_individual_points:
            sns.stripplot(
                y=plot_data, 
                color=point_color, 
                s=point_size,         # 使用 's' 而不是 'size'
                alpha=point_alpha,
                jitter=0.1,           # 稍微增加一点抖动，避免点完全在一条直线上
                ax=ax
            )

        # 设置图表标题和标签
        ax.set_title(plot_title, fontsize=16)
        ax.set_ylabel(y_axis_label, fontsize=12)
        ax.set_xlabel(column_name, fontsize=12)

        # 3. (新增) 设置自定义的Y轴范围
        if set_y_axis_range:
            print(f"设置自定义Y轴范围: ({y_axis_min}, {y_axis_max})")
            ax.set_ylim(y_axis_min, y_axis_max)

        # 调整布局
        plt.tight_layout()

        # 保存图表
        plt.savefig(output_file, dpi=300)
        
        print(f"图表已成功保存为: {output_file}")

    except FileNotFoundError:
        print(f"错误: 文件 '{input_file}' 未找到。")
    except Exception as e:
        print(f"生成图表时发生未知错误: {e}")

if __name__ == '__main__':
    create_violin_plot_v2()