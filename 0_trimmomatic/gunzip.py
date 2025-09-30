import os
import sys
import time

# --- 脚本配置 ---
# 1. 请在此处设置您想处理的文件后缀名
TARGET_EXTENSION = ".gz"

# 2. 请设置目标文件夹路径，使用 "." 代表当前脚本所在的文件夹
TARGET_DIRECTORY = "."
# --- 配置结束 ---


def find_files(directory: str, extension: str) -> list:
    """
    在指定目录中查找具有特定后缀名的文件。

    Args:
        directory (str): 要搜索的文件夹路径。
        extension (str): 目标文件的后缀名 (例如: ".gz")。

    Returns:
        list: 包含所有符合条件的文件名的列表。
    """
    print(f"[*] 正在 '{directory}' 文件夹中搜索 *{extension} 文件...")
    
    try:
        # 使用列表推导式，代码更简洁高效
        files = [f for f in os.listdir(directory) if f.endswith(extension)]
        return files
    except FileNotFoundError:
        print(f"[!] 错误：找不到文件夹 '{directory}'。请检查路径是否正确。")
        return []

def decompress_files(file_list: list, directory: str):
    """
    按顺序解压文件列表中的每一个文件。

    Args:
        file_list (list): 包含待解压文件名的列表。
        directory (str): 文件所在的目录路径。
    """
    total_files = len(file_list)
    if total_files == 0:
        print("[*] 没有找到需要处理的文件。")
        return

    print(f"[*] 共找到 {total_files} 个文件。准备开始解压...")
    
    start_time = time.time()

    # 使用 enumerate 来同时获取索引和文件名，方便打印进度
    for i, filename in enumerate(file_list):
        # 使用 os.path.join 来构建跨平台兼容的文件路径
        full_path = os.path.join(directory, filename)
        
        # 打印进度，\r 让光标回到行首，实现单行刷新效果
        # sys.stdout.flush() 用于确保立即输出
        progress_message = f"--> 正在处理: [{i + 1}/{total_files}] {filename}"
        sys.stdout.write('\r' + ' ' * 80) # 清除旧行
        sys.stdout.write('\r' + progress_message)
        sys.stdout.flush()

        # 执行解压命令
        os.system(f"gunzip {full_path}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    
    # 全部处理完后换行，避免覆盖最后一条进度
    print("\n") 
    print("=" * 30)
    print(f"[✓] 所有文件处理完毕！")
    print(f"[*] 总计用时: {elapsed_time:.2f} 秒。")
    print("=" * 30)


def main():
    """
    主执行函数
    """
    files_to_process = find_files(TARGET_DIRECTORY, TARGET_EXTENSION)
    decompress_files(files_to_process, TARGET_DIRECTORY)


# 推荐的Python脚本入口点
# 当这个.py文件被直接运行时，__name__ 的值是 '__main__'
# 如果它被其他脚本作为模块导入，则 __name__ 的值是模块名
if __name__ == "__main__":
    main()