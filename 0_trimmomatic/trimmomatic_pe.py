import os
import subprocess
import sys

# ==============================================================================
# ### CONFIGURATION (配置区域) ###
# 请在此处修改所有参数，无需通过命令行传参
# ==============================================================================

# 1. 路径设置
# ------------------------------------------------------------------------------
# 待处理文件所在的目录，默认为脚本当前所在目录 (os.getcwd())
# 如果需要指定其他目录，请修改为绝对路径，例如: "/home/user/data/raw_reads"
INPUT_DIRECTORY = os.getcwd()

# Trimmomatic 程序的 JAR 包路径
# 请确保该路径正确，如果不在当前目录，建议写绝对路径
TRIMMOMATIC_JAR_PATH = "trimmomatic-0.40.jar"

# 2. 文件后缀匹配设置
# ------------------------------------------------------------------------------
# 输入文件的后缀 (用于识别双端测序数据)
# 例如: 样本名为 "SampleA"，文件为 "SampleA_1.fq.gz"，则后缀为 "_1.fq.gz"
INPUT_FORWARD_SUFFIX = "1.fq.gz"  # 正向序列后缀 (R1)
INPUT_REVERSE_SUFFIX = "2.fq.gz"  # 反向序列后缀 (R2)

# Trimmomatic 输出的中间文件后缀 (通常不需要修改)
OUTPUT_PAIRED_FWD_SUFFIX = "_1_paired.fq"
OUTPUT_UNPAIRED_FWD_SUFFIX = "_1_unpaired.fq"
OUTPUT_PAIRED_REV_SUFFIX = "_2_paired.fq"
OUTPUT_UNPAIRED_REV_SUFFIX = "_2_unpaired.fq"

# 最终重命名后的输出文件后缀
FINAL_OUTPUT_FWD_SUFFIX = "_1.fq"
FINAL_OUTPUT_REV_SUFFIX = "_2.fq"

# 3. Trimmomatic 运行参数
# ------------------------------------------------------------------------------
# 线程数
THREADS = "30"

# 质量体系: "-phred33" 或 "-phred64"
PHRED_ENCODING = "-phred33"

# 具体的修剪步骤参数 (根据您的具体质控需求修改)
# 格式说明:
# LEADING/TRAILING: 去除首尾低质量碱基
# SLIDINGWINDOW: 滑动窗口修剪
# HEADCROP: 硬切除头部碱基
# MINLEN: 最小保留长度
TRIM_STEPS = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36"


# ==============================================================================
# ### END CONFIGURATION (配置结束) ###
# 下方代码为核心逻辑，非专业人士请勿修改
# ==============================================================================

def find_paired_reads(directory, f_suffix, r_suffix):
    """
    根据指定的后缀查找配对的测序文件。

    Args:
        directory (str): 搜索文件的目录路径。
        f_suffix (str): 正向序列文件的后缀 (R1)。
        r_suffix (str): 反向序列文件的后缀 (R2)。

    Returns:
        list: 包含配对样本基础名称 (basename) 的排序列表。
              如果未找到配对，返回空列表。
    """
    base_names = set()
    
    # 获取目录下所有文件列表
    try:
        all_files = os.listdir(directory)
    except FileNotFoundError:
        print(f"错误: 找不到目录 '{directory}'")
        return []

    print(f"正在目录 '{directory}' 中搜索后缀为 '{f_suffix}' 和 '{r_suffix}' 的文件...")

    # 将文件列表转换为集合以提高查找效率
    file_set = set(all_files)

    for filename in all_files:
        # 检查是否为正向文件
        if filename.endswith(f_suffix):
            # 获取样本基础名称 (去除后缀)
            base_name = filename[:-len(f_suffix)]

            # 构建预期的反向文件名
            expected_r_file = base_name + r_suffix

            # 检查对应的反向文件是否存在
            if expected_r_file in file_set:
                base_names.add(base_name)
            else:
                print(f"  警告: 发现正向文件 '{filename}' 但缺失对应的反向文件 '{expected_r_file}'。跳过。")

    if not base_names:
        print("  未找到完整的配对文件。")
    else:
        print(f"  共找到 {len(base_names)} 对潜在的配对样本。")

    return sorted(list(base_names))


def run_trimmomatic_workflow():
    """
    主执行函数：查找文件 -> 运行 Trimmomatic -> 清理中间文件 -> 重命名最终结果。
    """
    
    # 1. 获取配对样本
    sample_base_names = find_paired_reads(INPUT_DIRECTORY, INPUT_FORWARD_SUFFIX, INPUT_REVERSE_SUFFIX)

    if not sample_base_names:
        print("\n未找到符合条件的双端测序文件，程序退出。")
        sys.exit(1)

    print(f"\n即将处理以下 {len(sample_base_names)} 个样本: {sample_base_names}")

    # 记录处理失败的样本
    failed_samples = []

    # 2. 循环运行 Trimmomatic
    print("\n--- 开始运行 Trimmomatic 质控流程 ---")
    
    for base in sample_base_names:
        # 构建输入文件路径
        input_f_path = os.path.join(INPUT_DIRECTORY, base + INPUT_FORWARD_SUFFIX)
        input_r_path = os.path.join(INPUT_DIRECTORY, base + INPUT_REVERSE_SUFFIX)

        # 构建 Trimmomatic 输出路径
        out_pair_f = os.path.join(INPUT_DIRECTORY, base + OUTPUT_PAIRED_FWD_SUFFIX)
        out_unpair_f = os.path.join(INPUT_DIRECTORY, base + OUTPUT_UNPAIRED_FWD_SUFFIX)
        out_pair_r = os.path.join(INPUT_DIRECTORY, base + OUTPUT_PAIRED_REV_SUFFIX)
        out_unpair_r = os.path.join(INPUT_DIRECTORY, base + OUTPUT_UNPAIRED_REV_SUFFIX)

        # 构建命令列表
        cmd = [
            "java", "-jar", TRIMMOMATIC_JAR_PATH, "PE",
            "-threads", THREADS,
            PHRED_ENCODING,
            input_f_path, input_r_path,
            out_pair_f, out_unpair_f,
            out_pair_r, out_unpair_r
        ]
        
        # 添加修剪参数 (使用 split 将字符串转为列表参数)
        cmd.extend(TRIM_STEPS.split())

        print(f"\n正在处理样本: {base}")
        # 打印命令供调试 (将列表连接为字符串显示)
        print(f"  执行指令: {' '.join(cmd)}")

        try:
            # Python 3.6 兼容写法:
            # use stdout=subprocess.PIPE, stderr=subprocess.PIPE 代替 capture_output=True (3.7+)
            # use universal_newlines=True 代替 text=True (3.7+)
            process = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
                universal_newlines=True, 
                cwd=INPUT_DIRECTORY
            )
            
            # 打印运行日志
            if process.stdout:
                print(f"  [Trimmomatic STDOUT]:\n{process.stdout.strip()}")
            if process.stderr:
                # Trimmomatic 通常将统计信息输出到 stderr
                print(f"  [Trimmomatic STDERR]:\n{process.stderr.strip()}")
            
            print(f"  样本 {base} 处理成功。")

        except subprocess.CalledProcessError as e:
            print(f"  错误: 样本 '{base}' 处理失败!")
            print(f"  退出代码: {e.returncode}")
            print(f"  错误信息: {e.stderr.strip()}")
            failed_samples.append(base)
        except OSError as e:
            # 捕获找不到 java 命令或文件路径错误的情况
            print(f"  系统错误: 无法执行命令。请检查是否安装了 Java 或 JAR 路径是否正确。")
            print(f"  详细信息: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"  发生未知错误: {e}")
            failed_samples.append(base)

    # 3. 清理与重命名
    print("\n--- 开始清理中间文件并重命名 ---")

    for base in sample_base_names:
        if base in failed_samples:
            print(f"\n跳过失败样本的清理步骤: {base}")
            continue

        print(f"\n正在整理样本: {base}")

        # 定义需要删除的非配对文件路径
        files_to_remove = [
            os.path.join(INPUT_DIRECTORY, base + OUTPUT_UNPAIRED_FWD_SUFFIX),
            os.path.join(INPUT_DIRECTORY, base + OUTPUT_UNPAIRED_REV_SUFFIX)
        ]

        # 定义重命名规则: (旧文件名, 新文件名)
        rename_map = [
            (
                os.path.join(INPUT_DIRECTORY, base + OUTPUT_PAIRED_FWD_SUFFIX),
                os.path.join(INPUT_DIRECTORY, base + FINAL_OUTPUT_FWD_SUFFIX)
            ),
            (
                os.path.join(INPUT_DIRECTORY, base + OUTPUT_PAIRED_REV_SUFFIX),
                os.path.join(INPUT_DIRECTORY, base + FINAL_OUTPUT_REV_SUFFIX)
            )
        ]

        # 执行删除
        for f_path in files_to_remove:
            try:
                if os.path.exists(f_path):
                    os.remove(f_path)
                    print(f"  已删除中间文件: {os.path.basename(f_path)}")
            except OSError as e:
                print(f"  删除失败 '{os.path.basename(f_path)}': {e}")

        # 执行重命名
        for old_path, new_path in rename_map:
            try:
                if os.path.exists(old_path):
                    os.rename(old_path, new_path)
                    print(f"  已重命名: {os.path.basename(old_path)} -> {os.path.basename(new_path)}")
                else:
                    print(f"  警告: 找不到源文件 {os.path.basename(old_path)} (可能是质控后数据量为0?)")
            except OSError as e:
                print(f"  重命名失败 '{os.path.basename(old_path)}': {e}")

    print("\n--- 所有任务已完成 ---")

if __name__ == "__main__":
    run_trimmomatic_workflow()