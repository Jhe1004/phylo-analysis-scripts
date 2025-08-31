import os
import re
from Bio import SeqIO

def clean_name(name):
    """
    清理名称字符串，将空格和特殊字符替换为下划线。
    """
    # 定义需要被替换的特殊字符
    # 使用正则表达式，'[' 和 ']' 需要转义
    special_chars = r" |\(|\)|\[|\]"
    
    # 第一步: 使用正则表达式将所有特殊字符替换为下划线
    cleaned = re.sub(special_chars, "_", name)
    
    # 第二步: 将连续的多个下划线替换为单个下划线
    cleaned = re.sub(r"__+", "_", cleaned)
    
    # 第三步: 移除可能出现在开头或结尾的下划线
    cleaned = cleaned.strip('_')
    
    return cleaned

def extract_and_concatenate_cds():
    """
    主函数，执行CDS提取和合并的整个流程。
    """
    # --- 1. 设置输入和输出文件夹 ---
    input_directory = "input_gb_files"
    output_directory = "output_cds_fasta"

    # --- 2. 检查并创建文件夹 ---
    if not os.path.isdir(input_directory):
        print(f"错误：输入文件夹 '{input_directory}' 不存在。")
        print("请先创建该文件夹，并将您的.gb文件放入其中。")
        return

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print(f"已创建输出文件夹：'{output_directory}'")

    # --- 3. 遍历输入文件夹中的所有GenBank文件 ---
    print("\n开始处理GenBank文件...")
    processed_files = 0
    for filename in os.listdir(input_directory):
        if filename.endswith(".gb") or filename.endswith(".gbk"):
            input_filepath = os.path.join(input_directory, filename)
            
            # 从原始文件名中获取基础名（不含扩展名）
            original_base_name = os.path.splitext(filename)[0]
            
            # **【核心修改】**
            # 使用我们新的clean_name函数来处理名称
            species_name = clean_name(original_base_name)
            
            print(f"  -> 正在处理: '{original_base_name}' -> 清理后名称: '{species_name}'")

            concatenated_cds = ""
            
            try:
                for record in SeqIO.parse(input_filepath, "genbank"):
                    for feature in record.features:
                        if feature.type == "CDS":
                            if 'translation' in feature.qualifiers and 'pseudo' not in feature.qualifiers:
                                sequence = str(feature.location.extract(record).seq)
                                if len(sequence) % 3 == 0:
                                    concatenated_cds += sequence
                                else:
                                    print(f"    - 警告: 物种 {species_name} 中一个CDS序列长度不是3的倍数，已跳过。")
            except Exception as e:
                print(f"    - 错误: 处理文件 {filename} 时发生错误: {e}")
                continue

            # --- 4. 将合并后的序列写入FASTA文件 ---
            if concatenated_cds:
                # 使用清理后的名称来命名输出文件
                output_filename = f"{species_name}.fasta"
                output_filepath = os.path.join(output_directory, output_filename)
                
                with open(output_filepath, "w") as f_out:
                    # **【核心修改】**
                    # FASTA文件的头部也使用清理后的名称
                    f_out.write(f">{species_name}\n")
                    
                    for i in range(0, len(concatenated_cds), 60):
                        f_out.write(concatenated_cds[i:i+60] + "\n")
                
                print(f"    - 成功: 已将合并的CDS序列保存到 {output_filepath}")
                processed_files += 1
            else:
                print(f"    - 警告: 在物种 {species_name} 中未找到任何有效的CDS序列。")

    print(f"\n处理完成！总共成功处理并生成了 {processed_files} 个物种的FASTA文件。")


if __name__ == "__main__":
    extract_and_concatenate_cds()