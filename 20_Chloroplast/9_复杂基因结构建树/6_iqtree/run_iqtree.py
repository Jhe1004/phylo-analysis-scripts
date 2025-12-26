import subprocess
import sys
import os

# --- 您需要修改的设置 ---

# 1. 定义序列比对文件的名字
# !!! 把这里的文件名改成你自己的 FASTA 文件名 !!!
alignment_file = "my_alignment.fasta"

# 2. 定义分区文件的名字
# 确保这个名字和你的分区文件名一致
partition_file = "partitions.txt"

# 3. 定义iqtree命令参数
# 您可以根据需要修改这些参数
# -m MFP: 自动为每个分区寻找最佳模型 (ModelFinder Plus)
# -bb 1000: 运行1000次 Ultrafast Bootstrap
# -nt AUTO: 自动检测可用的CPU核心数
IQTRRE_ARGS = [
    "-m", "MFP",
    "-bb", "1000",
    "-nt", "AUTO"
]

# --- 脚本主程序 (下面一般不需要修改) ---

def run_iqtree():
    # 1. 检查所需文件是否存在
    if not os.path.exists(alignment_file):
        print(f"错误: 找不到比对文件: {alignment_file}")
        print("请确保您在脚本顶部的 'alignment_file' 变量中设置了正确的文件名。")
        sys.exit(1)
        
    if not os.path.exists(partition_file):
        print(f"错误: 找不到分区文件: {partition_file}")
        print(f"请确保一个名为 '{partition_file}' 的有效分区文件在此目录下。")
        print("根据您的示例，它看起来像这样：")
        print("DNA, part1 = 1-2249")
        print("DNA, part2 = 2250-2848")
        sys.exit(1)
    
    # 2. 构建完整的iqtree命令
    # 我们使用 -spp 参数来指定分区文件
    iqtree_command = [
        "iqtree",
        "-s", alignment_file,
        "-spp", partition_file
    ]
    
    # 把您在上面定义的其他参数加进来
    iqtree_command.extend(IQTRRE_ARGS)
    
    print("--- 准备运行 IQ-TREE ---")
    print("命令: " + " ".join(iqtree_command))
    print("--------------------------")
    
    # 3. 运行命令
    try:
        # 运行命令，并等待它完成
        subprocess.run(iqtree_command, check=True)
        print("\n--- IQ-TREE 运行成功! ---")
        
    except subprocess.CalledProcessError as e:
        # 如果iqtree返回了错误码
        print(f"\n--- IQ-TREE 运行失败 ---")
        print(f"错误信息: {e}")
        
    except FileNotFoundError:
        # 如果系统中找不到 'iqtree' 命令
        print("\n--- 错误: 'iqtree' 命令未找到 ---")
        print("请确保您已经正确安装了 IQ-TREE，")
        print("并且 'iqtree' 可执行文件在您的系统 PATH 环境变量中。")

if __name__ == "__main__":
    run_iqtree()