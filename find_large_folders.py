import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

def get_directory_stats(directory):
    """
    递归计算目录的总大小和总文件数。
    返回一个元组 (total_size_bytes, total_file_count)。
    """
    total_size = 0
    total_count = 0
    try:
        for entry in os.scandir(directory):
            if entry.is_file(follow_symlinks=False):
                total_size += entry.stat(follow_symlinks=False).st_size
                total_count += 1
            elif entry.is_dir(follow_symlinks=False):
                dir_size, dir_count = get_directory_stats(entry.path)
                total_size += dir_size
                total_count += dir_count
    except (PermissionError, FileNotFoundError):
        pass
    return total_size, total_count

def find_problematic_directories(root_path, min_size_mb, min_file_count):
    """
    使用多线程查找指定路径下体积过大或文件过多的文件夹。
    """
    min_size_bytes = min_size_mb * 1024 * 1024
    problematic_directories = []
    
    directories_to_scan = [dirpath for dirpath, _, _ in os.walk(root_path)]

    print(f"准备分析 {len(directories_to_scan)} 个文件夹...")

    with ThreadPoolExecutor() as executor:
        future_to_dir = {executor.submit(get_directory_stats, d): d for d in directories_to_scan}

        for future in as_completed(future_to_dir):
            directory = future_to_dir[future]
            try:
                size, count = future.result()
                reasons = []
                if size >= min_size_bytes:
                    reasons.append("体积过大")
                if count >= min_file_count:
                    reasons.append("文件过多")

                if reasons:
                    problematic_directories.append({
                        "path": directory,
                        "size": size,
                        "count": count,
                        "reason": " 和 ".join(reasons)
                    })
            except Exception as exc:
                print(f"'{directory}' 生成了一个异常: {exc}", file=sys.stderr)

    return problematic_directories

def filter_nested_results(directories):
    """
    过滤掉结果中的父文件夹，只保留最深层的“问题”文件夹。
    """
    if not directories:
        return []
    
    paths = [os.path.normpath(d['path']) for d in directories]
    paths_to_remove = set()

    for p1 in paths:
        for p2 in paths:
            if p1 == p2:
                continue
            if p2.startswith(p1 + os.sep):
                paths_to_remove.add(p1)
                break
    
    final_results = [d for d in directories if os.path.normpath(d['path']) not in paths_to_remove]
    return final_results

def format_size(size_bytes):
    """
    将字节大小格式化为更易读的单位 (B, KB, MB, GB)。
    """
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024**2:
        return f"{size_bytes/1024:.2f} KB"
    elif size_bytes < 1024**3:
        return f"{size_bytes/1024**2:.2f} MB"
    else:
        return f"{size_bytes/1024**3:.2f} GB"

if __name__ == "__main__":
    # --- 参数配置 ---
    target_directory = '.'
    size_threshold_mb = 20
    file_count_threshold = 1000

    print(f"正在扫描 '{os.path.abspath(target_directory)}' 及其所有子文件夹...")
    print(f"筛选标准: 体积超过 {size_threshold_mb} MB 或 文件数超过 {file_count_threshold} 个\n")

    try:
        all_problem_dirs = find_problematic_directories(target_directory, size_threshold_mb, file_count_threshold)

        print(f"\n初步找到 {len(all_problem_dirs)} 个有问题文件夹，正在过滤冗余结果...")
        final_dirs = filter_nested_results(all_problem_dirs)
        print("-" * 40)

        # --- 分类和排序 ---
        large_size_dirs = []
        many_files_dirs = []

        for d in final_dirs:
            if "体积过大" in d['reason']:
                large_size_dirs.append(d)
            if "文件过多" in d['reason']:
                many_files_dirs.append(d)
        
        # 按大小降序排序
        large_size_dirs.sort(key=lambda item: item['size'], reverse=True)
        # 按文件数降序排序
        many_files_dirs.sort(key=lambda item: item['count'], reverse=True)

        # --- 分类打印结果 ---
        if not large_size_dirs and not many_files_dirs:
            print("\n未找到需要关注的文件夹。")
        else:
            # 打印体积过大的文件夹
            if large_size_dirs:
                print(f"\n--- 发现 {len(large_size_dirs)} 个【体积过大】的文件夹 (按大小排序) ---")
                for info in large_size_dirs:
                    path = info['path']
                    size_str = format_size(info['size'])
                    count_str = f"{info['count']} 个文件"
                    print(f"  - 路径: {path}\n    | 大小: {size_str}, 文件数: {count_str}\n")
            
            # 打印文件过多的文件夹
            if many_files_dirs:
                print(f"\n--- 发现 {len(many_files_dirs)} 个【文件过多】的文件夹 (按文件数排序) ---")
                for info in many_files_dirs:
                    path = info['path']
                    size_str = format_size(info['size'])
                    count_str = f"{info['count']} 个文件"
                    print(f"  - 路径: {path}\n    | 大小: {size_str}, 文件数: {count_str}\n")

    except FileNotFoundError:
        print(f"错误: 找不到指定的目录 '{target_directory}'", file=sys.stderr)
    except Exception as e:
        print(f"发生未知错误: {e}", file=sys.stderr)