import sys
import os
from collections import OrderedDict

def parse_fasta(filename):
    """
    解析FASTA文件，返回一个包含(name, sequence)元组的列表以保持顺序。
    """
    sequences = []
    current_seq = []
    current_name = None
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_name:
                        sequences.append((current_name, "".join(current_seq)))
                    current_name = line[1:].strip()
                    current_seq = []
                else:
                    # 移除序列中可能存在的空格
                    current_seq.append(line.replace(" ", ""))
        
        # 添加最后一个序列
        if current_name:
            sequences.append((current_name, "".join(current_seq)))
            
    except FileNotFoundError:
        print(f"错误：找不到文件 {filename}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"读取FASTA文件时出错: {e}", file=sys.stderr)
        return None
        
    if not sequences:
        print(f"警告：文件 {filename} 中未找到序列。", file=sys.stderr)
        return None
        
    return sequences

def format_seq_chunk(seq_chunk, group_size=10):
    """
    将序列块格式化为每 'group_size' 个字符加一个空格。
    """
    parts = []
    for i in range(0, len(seq_chunk), group_size):
        parts.append(seq_chunk[i:i+group_size])
    return " ".join(parts)

def convert_fasta_to_phylip(fasta_file, phylip_file, is_batch=False):
    """
    将FASTA格式的比对文件转换为严格的、与示例文件格式一致的PHYLIP文件。
    is_batch: True时，减少打印输出，用于批量模式。
    """
    if not is_batch:
        print(f"--- 模式：单文件转换 ---")
        print(f"开始转换 {fasta_file}...")
    
    sequences = parse_fasta(fasta_file)
    
    if not sequences:
        # 在批量模式下，错误会由 run_batch_conversion 捕获
        raise Exception("未能解析FASTA文件或文件为空。")

    # --- 1. 验证数据 ---
    num_taxa = len(sequences)
    if num_taxa == 0:
        raise Exception("FASTA文件中未找到序列。")
        
    alignment_length = len(sequences[0][1])
    # 检查所有序列是否等长
    for name, seq in sequences[1:]:
        if len(seq) != alignment_length:
            error_msg = f"错误：序列长度不一致。序列 {sequences[0][0]} ({alignment_length}bp) 与 序列 {name} ({len(seq)}bp) 长度不同。"
            if not is_batch:
                print(error_msg, file=sys.stderr)
                print("请确保输入的是一个比对(alignment)文件。", file=sys.stderr)
            raise Exception(error_msg)

    # --- 2. 确定格式参数 ---
    max_name_len = 0
    for name, _ in sequences:
        if len(name) > max_name_len:
            max_name_len = len(name)
            
    # 名称列的宽度 = 最长名称 + 2个空格
    NAME_COLUMN_WIDTH = max_name_len + 2
    CHUNK_SIZE = 60
    GROUP_SIZE = 10
    
    if not is_batch:
        print(f"检测到 {num_taxa} 个序列，比对长度 {alignment_length} bp。")
        print(f"最长名称为 {max_name_len} 字符，名称列宽将设置为 {NAME_COLUMN_WIDTH}。")

    # --- 3. 写入PHYLIP文件 ---
    try:
        with open(phylip_file, 'w', encoding='utf-8') as f:
            # 写入PHYLIP头部
            f.write(f" {num_taxa}   {alignment_length}\n")
            
            # 循环遍历序列的每个区块
            for i in range(0, alignment_length, CHUNK_SIZE):
                if i > 0:
                    f.write("\n")  # 在区块之间添加一个空行
                
                start_index = i
                end_index = min(i + CHUNK_SIZE, alignment_length)
                
                # 遍历每个序列，写入当前区块
                for name, seq in sequences:
                    seq_chunk = seq[start_index:end_index]
                    formatted_chunk = format_seq_chunk(seq_chunk, GROUP_SIZE)
                    
                    if i == 0:
                        # 第一个区块：写入名称 + 序列
                        padded_name = name.ljust(NAME_COLUMN_WIDTH)
                        f.write(f"{padded_name}{formatted_chunk}\n")
                    else:
                        # 后续区块：写入空格 + 序列
                        padding_spaces = " " * NAME_COLUMN_WIDTH
                        f.write(f"{padding_spaces}{formatted_chunk}\n")
                        
        if not is_batch:
            print(f"\n成功！文件已保存为 {phylip_file}")
        else:
            # 在批量模式下，简化成功信息
            print(f"  [成功] {os.path.basename(phylip_file)} 已保存。")

    except Exception as e:
        # 将异常抛出，由调用者处理
        raise Exception(f"写入PHYLIP文件时出错: {e}")


def run_batch_conversion():
    """
    在当前文件夹中查找所有FASTA文件并批量转换为PHYLIP格式。
    """
    print(f"--- 模式：批量转换 (当前文件夹) ---")
    
    FASTA_EXTENSIONS = ('.fasta', '.fa', '.fas', '.fna')
    current_directory = os.getcwd()
    files_processed = 0
    files_skipped = 0
    files_failed = 0

    print(f"正在扫描: {current_directory}\n")

    for filename in os.listdir(current_directory):
        # 检查文件扩展名
        if filename.lower().endswith(FASTA_EXTENSIONS):
            input_file = os.path.join(current_directory, filename)
            # 输出文件名 (例如: my_align.fasta -> my_align.phy)
            output_file = os.path.join(current_directory, os.path.splitext(filename)[0] + ".phy")
            
            # 检查输出文件是否已存在
            if os.path.exists(output_file):
                print(f"[跳过] {filename} (输出文件 {os.path.basename(output_file)} 已存在)")
                files_skipped += 1
                continue
                
            print(f"[处理中] {filename} -> {os.path.basename(output_file)}")
            
            try:
                # 调用核心转换函数
                convert_fasta_to_phylip(input_file, output_file, is_batch=True)
                files_processed += 1
            except Exception as e:
                print(f"  [错误] 处理 {filename} 时失败: {e}", file=sys.stderr)
                files_failed += 1
            
            print("-" * 20) # 添加分隔符

    print(f"\n--- 批量处理完成 ---")
    print(f"成功转换: {files_processed} 个文件")
    print(f"跳过: {files_skipped} 个文件 (已存在)")
    print(f"失败: {files_failed} 个文件")
    
    if (files_processed + files_skipped + files_failed) == 0:
        print("未在当前文件夹中找到任何 FASTA 文件。")


# --- 脚本执行入口 ---
if __name__ == "__main__":
    
    # 模式1: 批量处理 (无参数)
    # 用法: python fasta_to_phylip.py
    if len(sys.argv) == 1:
        run_batch_conversion()
        
    # 模式2: 单文件转换 (2个参数)
    # 用法: python fasta_to_phylip.py <input.fasta> <output.phy>
    elif len(sys.argv) == 3:
        input_fasta = sys.argv[1]
        output_phylip = sys.argv[2]
        
        if not os.path.exists(input_fasta):
            print(f"错误：输入文件不存在 {input_fasta}", file=sys.stderr)
            sys.exit(1)
            
        try:
            convert_fasta_to_phylip(input_fasta, output_phylip, is_batch=False)
        except Exception as e:
            print(f"\n转换时发生错误: {e}", file=sys.stderr)
            sys.exit(1)
            
    # 模式3: 用法错误
    else:
        print("用法错误。", file=sys.stderr)
        print("\n请选择一种模式:", file=sys.stderr)
        print("  1. 批量模式: python fasta_to_phylip.py", file=sys.stderr)
        print("     (自动转换当前文件夹中所有的 .fasta, .fa, .fas, .fna 文件)", file=sys.stderr)
        print("\n  2. 单文件模式: python fasta_to_phylip.py <input.fasta> <output.phy>", file=sys.stderr)
        sys.exit(1)