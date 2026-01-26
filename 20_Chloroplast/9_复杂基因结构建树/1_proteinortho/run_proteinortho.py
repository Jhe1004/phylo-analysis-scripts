import os
import glob
import sys  # 导入 sys 库用于退出脚本

# --- 脚本配置 ---
# 定义您要搜索的FASTA文件的扩展名
# 这应该匹配您上一个脚本 (split_genome_by_feature.py) 生成的文件
FASTA_EXTENSION = ".fasta"
# ---

print("开始运行 proteinortho (DNA 模式) 脚本...")

# 1. 获取当前目录下所有 .fasta 文件
#    (注意：这里假设您每个物种的片段在一个单独的FASTA文件中)
fasta_files = glob.glob(f"*{FASTA_EXTENSION}")

# 检查是否找到了文件
if not fasta_files:
    print(f"错误：在当前目录未找到 *{FASTA_EXTENSION} 文件。")
    print("请确保您上一步生成的 FASTA 文件（每个物种一个）在此文件夹中。")
    sys.exit(1)  # 如果没找到文件，退出脚本

print(f"找到 {len(fasta_files)} 个 *{FASTA_EXTENSION} 文件准备进行比较。")

# 2. 清理步骤 (已修改)
#    您之前的脚本是为蛋白质 (.pep) 准备的，用于移除 '*' (终止密码子)。
#    对于DNA (.fasta) 文件，这个清理步骤不是必需的，因为DNA序列没有 '*'。
#    因此，我们跳过这个循环。如果您需要对DNA序列进行其他清理（例如移除'N'之外的特殊字符），
#    可以在这里添加。

# print("开始清理文件 (此步骤对DNA文件通常不是必需的)...")
# for fasta_file in fasta_files:
#     temp_file = fasta_file + '.tmp'
#     
#     with open(fasta_file, 'r') as f_in, open(temp_file, 'w') as f_out:
#         for line in f_in:
#             # --- 这是您之前的逻辑，针对蛋白质，已注释掉 ---
#             # if not line.startswith('>'):
#             #     line = line.replace('*', '')
#             f_out.write(line)
#     
#     os.replace(temp_file, fasta_file)
#     print(f"已“清理”文件: {fasta_file}")
# print("所有文件清理完毕。")


# 3. 运行 Proteinortho (已修改)
print("\n即将运行 proteinortho (使用 blastn DNA模式)...")

# 构建命令
# -p=blastn 告诉 proteinortho 使用 blastn (DNA vs DNA)
# f"*{FASTA_EXTENSION}" 会自动匹配所有 .fasta 文件
command = f"proteinortho -p=blastn *{FASTA_EXTENSION}"

print(f"执行命令: {command}")

# 运行命令
# 注意：请确保 'proteinortho' 和 'blastn' 已经安装
# 并且在您系统的 PATH 环境变量中。
try:
    os.system(command)
    print("\nproteinortho 运行完成。")
    print("请检查 'myproject.proteinortho' 和 'myproject.proteinortho.blast-graph' 等输出文件。")
except Exception as e:
    print(f"\n*** 运行 proteinortho 时出错: {e} ***")
    print("请确保 proteinortho 和 BLAST (blastn) 已正确安装并配置在系统路径中。")