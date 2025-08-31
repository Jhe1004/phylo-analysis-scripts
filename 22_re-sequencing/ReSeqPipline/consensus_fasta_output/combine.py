import os
import glob

def combine_fasta_files(output_filename="combined_sequences.fasta"):
    """
    将当前目录中的所有 .fasta 文件合并成一个单一的 FASTA 文件。
    合并后的文件中，每条序列的名称来源于其原始文件名。

    参数:
    output_filename (str): 合并后输出的FASTA文件名。
    """
    
    # 获取当前目录下所有 .fasta 文件
    # 使用 os.path.abspath 获取绝对路径以避免潜在问题
    # 使用 os.getcwd() 获取当前工作目录
    current_directory = os.getcwd()
    fasta_files_pattern = os.path.join(current_directory, "*.fasta")
    fasta_files = glob.glob(fasta_files_pattern)

    # 构建输出文件的绝对路径
    absolute_output_filename = os.path.join(current_directory, output_filename)

    # 从文件列表中排除输出文件自身，以防脚本重复运行时读取自己
    # glob.glob返回的是包含路径的文件名列表，所以比较时也要用完整路径
    if absolute_output_filename in fasta_files:
        fasta_files.remove(absolute_output_filename)
    
    if not fasta_files:
        print("在当前目录中没有找到 .fasta 文件。")
        return

    try:
        with open(absolute_output_filename, 'w') as outfile:
            for filepath in fasta_files:
                # 提取文件名（不含扩展名）作为序列头部
                # os.path.basename 获取文件名（含扩展名）
                # os.path.splitext 分割文件名和扩展名
                filename_without_extension = os.path.splitext(os.path.basename(filepath))[0]
                
                # 创建新的FASTA头部
                new_header = f">{filename_without_extension}\n"
                
                sequence_parts = []
                try:
                    with open(filepath, 'r') as infile:
                        for line in infile:
                            line = line.strip()  # 去除行首尾的空白字符
                            if not line:  # 跳过空行
                                continue
                            if line.startswith(">"):  # 跳过原始文件的FASTA头部行
                                continue
                            sequence_parts.append(line) # 添加序列部分
                    
                    if not sequence_parts:
                        print(f"警告: 文件 {os.path.basename(filepath)} 中没有找到有效的序列数据。已跳过。")
                        continue

                    full_sequence = "".join(sequence_parts)
                    
                    # 将新的头部和序列写入到合并文件中
                    outfile.write(new_header)
                    outfile.write(full_sequence + "\n") # 每条序列后加一个换行符
                
                except Exception as e_inner:
                    print(f"处理文件 {os.path.basename(filepath)} 时发生错误: {e_inner}")

        print(f"成功将所有 .fasta 文件合并到 {output_filename}")

    except Exception as e_outer:
        print(f"创建或写入输出文件 {output_filename} 时发生错误: {e_outer}")

if __name__ == "__main__":
    # 你可以将希望的输出文件名传递给函数，例如:
    # combine_fasta_files("my_custom_combined_name.fasta")
    combine_fasta_files()

    # --- 使用示例 ---
    # 1. 将此脚本保存为 .py 文件（例如: merge_fasta.py）到你的 .fasta 文件所在的文件夹。
    # 2. 打开终端或命令行。
    # 3. 导航到该文件夹 (例如: cd /path/to/your/fasta_files_folder)。
    # 4. 运行脚本: python merge_fasta.py
    # 5. 脚本执行完毕后，会在同一文件夹下生成一个名为 "combined_sequences.fasta" 的文件，
    #    其中包含了所有原始 .fasta 文件中的序列。