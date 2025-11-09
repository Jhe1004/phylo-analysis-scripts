import os
import glob
import sys
import multiprocessing

# --- 用户配置 ---

# 1. 查找包含FASTA文件的子文件夹的通配符
SUBDIR_PATTERN = "*_groups"

# 2. 目标文件的扩展名
FILE_TAG = ".fasta"

# 3. CPU 和线程管理
# 自动检测CPU核心数
try:
    TOTAL_CORES = multiprocessing.cpu_count()
except NotImplementedError:
    TOTAL_CORES = 4  # 默认值
    
# 并行运行的MAFFT进程数 (设置为您的总核心数)
NUM_CPU = TOTAL_CORES
# 每个MAFFT进程使用的线程数 (对于大量小文件, 1是最高效的)
THREAD_PER_JOB = 1

print(f"检测到 {TOTAL_CORES} 个CPU核心。")
print(f"将使用 {NUM_CPU} 个并行进程，每个进程 {THREAD_PER_JOB} 个线程。")


# --- 脚本正文 ---

def get_file_list_in_subdir():
    """
    获取当前子文件夹中符合目标扩展名, 且尚未被比对的文件列表。
    """
    file_name = []
    
    # 获取当前文件夹的名字, 用于排除最终的concat文件
    current_dir_name = os.path.basename(os.getcwd())
    output_concat_file = f"{current_dir_name}.fasta"

    for each in os.listdir(os.getcwd()):
        # 1. 必须是 .fasta 文件
        # 2. 不能是已经mafft过的
        # 3. 不能是最终连接好的
        if (each.endswith(FILE_TAG) and 
            not each.endswith("_maffted.fasta") and 
            each != output_concat_file):
            
            file_name.append(each)
            
    return file_name

def run_mafft(each_file):
    """
    对单个文件运行MAFFT（包含方向调整和sed清理）。
    """
    output_file = each_file.replace(FILE_TAG, "") + "_maffted.fasta"
    
    # 保留您命令中的 --adjustdirection 和 sed
    command = (
        f"mafft --adjustdirection --auto --thread {THREAD_PER_JOB} "
        f"\"{each_file}\" | sed 's/_R_//g' > \"{output_file}\""
    )
    
    try:
        # print(f"  正在处理: {each_file}") # 如果需要详细日志, 取消此行注释
        os.system(command)
        # print(f"  处理完成: {each_file}")
    except Exception as e:
        print(f"  处理文件 {each_file} 时出错: {e}")

def main():
    """
    主函数：遍历所有子目录并批量运行MAFFT。
    """
    main_dir = os.getcwd()
    print(f"主目录: {main_dir}")
    
    # 1. 查找所有匹配的子文件夹
    subdirectories = [d for d in glob.glob(SUBDIR_PATTERN) if os.path.isdir(d)]
    
    if not subdirectories:
        print(f"错误: 在 {main_dir} 中未找到与 '{SUBDIR_PATTERN}' 匹配的子文件夹。")
        sys.exit(1)
        
    print(f"找到 {len(subdirectories)} 个子文件夹: {', '.join(subdirectories)}\n")
    
    # 2. 遍历每个子文件夹
    for subdir in subdirectories:
        subdir_path = os.path.join(main_dir, subdir)
        print(f"--- 正在进入: {subdir} ---")
        
        # 切换到子目录
        try:
            os.chdir(subdir_path)
        except Exception as e:
            print(f"  错误: 无法进入文件夹 {subdir_path}: {e}")
            continue

        # 3. 查找需要比对的FASTA文件
        file_list = get_file_list_in_subdir()
        
        if not file_list:
            print("  没有找到需要比对的 .fasta 文件。跳至下一个文件夹。\n")
            os.chdir(main_dir) # 切换回主目录
            continue
            
        print(f"  找到 {len(file_list)} 个文件准备进行MAFFT比对...")
        
        # 4. 使用多进程并行运行MAFFT
        with multiprocessing.Pool(processes=NUM_CPU) as pool:
            # 使用 pool.map 等待所有任务完成
            pool.map(run_mafft, file_list)

        print(f"  {subdir} 中的所有MAFFT任务已完成。\n")
        
        # 5. 切换回主目录
        os.chdir(main_dir)

    print("--- 所有批量MAFFT任务已完成 ---")

# --- 运行脚本 ---
if __name__ == "__main__":
    
    # 检查 mafft 和 sed 是否已安装 (一个简单的检查)
    if os.system("mafft --version > /dev/null 2>&1") != 0:
         # 适用于 Windows 的备用检查
        if os.system("mafft --version > NUL 2>&1") != 0:
            print("错误: 未在您的系统路径中找到 'mafft'。")
            print("请先安装 MAFFT 并将其添加到您的 PATH 环境变量中。")
            sys.exit(1)
            
    if os.system("sed --version > /dev/null 2>&1") != 0:
        if os.system("sed --version > NUL 2>&1") != 0:
            print("警告: 未找到 'sed' 命令。")
            print("如果您在 Windows (非WSL) 环境, MAFFT 产生的 '_R_' 标记可能不会被移除。")

    main()