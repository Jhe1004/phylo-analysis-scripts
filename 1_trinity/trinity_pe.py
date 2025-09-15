# -*- coding: utf-8 -*-

import os
import sys

# ======================= 用户配置区 =======================
# 请在此处修改您希望处理的文件后缀名列表。
# 脚本会按照列表中的顺序查找匹配的文件。
# 例如: ['_1.fq.gz', '_2.fq.gz'] or ['.fq', '.fastq']
VALID_EXTENSIONS = ['.fq', '.fastq'] 
# ==========================================================


def find_paired_end_files(directory, extensions):
    """
    在指定目录中查找成对的测序文件 (例如 sample_1.fq 和 sample_2.fq)。

    Args:
        directory (str): 要搜索的目录。
        extensions (list): 有效的文件后缀名列表。

    Returns:
        dict: 一个字典，键是样本的基础名称，值是对应的文件后缀名。
              例如: {'sampleA': '.fq', 'sampleB': '.fastq'}
    """
    print("开始扫描文件...")
    sample_files = {}
    files_in_dir = os.listdir(directory)
    
    for ext in extensions:
        # 查找所有以 "_1" 和指定后缀结尾的文件
        for filename in files_in_dir:
            if filename.endswith("_1" + ext):
                # 构造对应的 "_2" 文件名并检查它是否存在
                right_file = filename.replace("_1" + ext, "_2" + ext)
                if right_file in files_in_dir:
                    # 提取样本的基础名称 (例如从 'sampleA_1.fq' 得到 'sampleA')
                    base_name = filename.replace("_1" + ext, "")
                    if base_name not in sample_files:
                        print(f"  [成功] 找到样本: {base_name} (后缀: {ext})")
                        sample_files[base_name] = ext
                    else:
                        # 如果之前已经找到了一个带有不同后缀的同名文件，给出提示
                        print(f"  [警告] 样本 '{base_name}' 存在多种后缀的文件，将使用已找到的 '{sample_files[base_name]}'")
    
    if not sample_files:
        print("\n错误: 在当前目录中未找到任何符合命名规范的成对文件。")
        print(f"请确保您的文件名以 '_1' 和 '_2' 结尾，且后缀为 {VALID_EXTENSIONS} 中的一种。")
        sys.exit(1)  # 如果找不到任何文件，则退出脚本
        
    return sample_files

def main():
    """
    主执行函数
    """
    now_dir = os.getcwd()
    print(f"当前工作目录: {now_dir}\n")
    
    # 1. 查找所有符合条件的成对文件
    samples_to_process = find_paired_end_files(now_dir, VALID_EXTENSIONS)
    
    print("\n--------------------------------------------------")
    print("准备开始处理以下样本:")
    for sample, ext in samples_to_process.items():
        print(f"- {sample} (使用后缀: {ext})")
    print("--------------------------------------------------\n")

    # 2. 为每个样本构建并执行Trinity命令
    for sample_name, extension in samples_to_process.items():
        left_file = os.path.join(now_dir, f"{sample_name}_1{extension}")
        right_file = os.path.join(now_dir, f"{sample_name}_2{extension}")
        output_dir = os.path.join(now_dir, f"{sample_name}_trinity")
        
        # 使用 f-string 构建命令，更清晰易读
        # 注意：文件和目录路径用引号括起来，以防路径中包含空格
        command = (
            "singularity exec -e trinityrnaseq.v2.15.1.simg Trinity "
            f"--seqType fq "
            f"--left \"{left_file}\" "
            f"--right \"{right_file}\" "
            f"--max_memory 40G "
            f"--CPU 40 "
            f"--output \"{output_dir}\" "
            f"--full_cleanup"
        )
        
        print(f">>> 正在准备处理样本: {sample_name}")
        print(f"将要执行的命令是:\n{command}\n")
        

        os.system(command)
        
        print(f"--- 样本 {sample_name} 的命令已生成 ---\n")

if __name__ == "__main__":
    main()