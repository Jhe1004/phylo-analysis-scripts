import os
import matplotlib.pyplot as plt
import numpy as np
from sklearn.mixture import GaussianMixture

# 指定包含 .txt 文件的文件夹路径
folder_path = './'  # 替换为你文件夹的实际路径

# 创建输出文件夹来保存图像
output_folder = 'output_images'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# 遍历文件夹中的所有 .txt 文件
for filename in os.listdir(folder_path):
    if filename.endswith('.pep'):  # 只处理 .txt 文件
        file_path = os.path.join(folder_path, filename)
        
        # 读取文件中的数据
        with open(file_path, 'r') as file:
            data = [float(line.strip()) for line in file if 5 <= float(line.strip()) <= 100]

        # 如果有数据，进行处理
        if data:
            data = np.array(data).reshape(-1, 1)  # 将数据转换为二维数组，适应GMM

            # 尝试用1到4个正态分布来拟合数据，选择最优模型
            lowest_bic = np.inf
            best_gmm = None
            n_components_range = range(1, 4)
            bics = []

            for n_components in n_components_range:
                gmm = GaussianMixture(n_components=n_components)
                gmm.fit(data)
                bic = gmm.bic(data)
                bics.append(bic)
                if bic < lowest_bic:
                    lowest_bic = bic
                    best_gmm = gmm

            # 绘制直方图，设置为浅灰色
            plt.figure(figsize=(5, 3))
            count, bins, ignored = plt.hist(data, bins=30, density=True, alpha=0.6, color='#f5f5f5', edgecolor='black')  # 设置color为浅灰色
            
            # 使用最佳模型绘制多个正态分布，曲线颜色统一为红色
            means = best_gmm.means_.flatten()
            covariances = np.sqrt(best_gmm.covariances_).flatten()
            weights = best_gmm.weights_.flatten()

            x = np.linspace(min(data), max(data), 1000).reshape(-1, 1)
            for mean, cov, weight in zip(means, covariances, weights):
                plt.plot(x, weight * (1/(cov * np.sqrt(2 * np.pi))) * np.exp(- (x - mean)**2 / (2 * cov**2)),
                         linewidth=2, color='#e4a19f')  # 将曲线颜色设置为红色并加粗

            # 调整坐标轴刻度字体大小
            plt.tick_params(axis='both', which='major', labelsize=11)

            # 保存图像
            output_image_path = os.path.join(output_folder, f'{filename[:-4]}_gmm.png')  # 去掉 .txt 后缀
            plt.savefig(output_image_path, transparent=True)  # 保存为图片
            plt.close()  # 关闭图表，避免内存占用过多
        else:
            print(f"No data in the range 10-100 found in {filename}.")
