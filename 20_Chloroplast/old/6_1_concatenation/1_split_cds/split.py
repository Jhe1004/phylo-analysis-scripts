#!/usr/bin/env python3
import os

def read_fasta(file_path):
    """
    读取给定fasta文件，返回一个列表，每个元素为 (header, sequence) 的元组。
    """
    sequences = []
    header = None
    seq_lines = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # 如果不是第一次遇到序列头，需要先存储前一个序列
                if header and seq_lines:
                    sequences.append((header, ''.join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        # 文件结束后将最后一个序列加入列表
        if header and seq_lines:
            sequences.append((header, ''.join(seq_lines)))
    return sequences

def split_sequence_into_frames(seq):
    """
    将一条核苷酸序列根据密码子位置（1/2/3位）拆分成三个子序列。
    假设序列为 "ATGCGT..." ，则:
    frame1 = A   C   ...
    frame2 = T   G   ...
    frame3 = G   T   ...
    """
    frame1 = seq[0::3]
    frame2 = seq[1::3]
    frame3 = seq[2::3]
    return frame1, frame2, frame3

def process_fasta_file(fasta_file, out_dir1, out_dir2, out_dir3):
    """
    处理单个fasta文件，对每条序列进行三分拆分，并将结果写入分别位于out_dir1、out_dir2、out_dir3下的文件中。
    """
    # 读取原始fasta文件内容
    sequences = read_fasta(fasta_file)
    
    # 构建输出文件名。例如原文件为 "example.fasta"
    base_name = os.path.basename(fasta_file)
    if base_name.lower().endswith('.fasta'):
        base_name = base_name[:-6]  # 去掉'.fasta'
    elif base_name.lower().endswith('.fa'):
        base_name = base_name[:-3]  # 有些文件可能是'.fa'
    
    # 输出文件路径
    out_file1 = os.path.join(out_dir1, f"{base_name}_frame1.fasta")
    out_file2 = os.path.join(out_dir2, f"{base_name}_frame2.fasta")
    out_file3 = os.path.join(out_dir3, f"{base_name}_frame3.fasta")

    # 打开输出文件进行写入
    with open(out_file1, 'w') as f1, open(out_file2, 'w') as f2, open(out_file3, 'w') as f3:
        for header, seq in sequences:
            frame1, frame2, frame3 = split_sequence_into_frames(seq)
            # 写入frame1序列
            f1.write(f">{header}\n{frame1}\n")
            # 写入frame2序列
            f2.write(f">{header}\n{frame2}\n")
            # 写入frame3序列
            f3.write(f">{header}\n{frame3}\n")

def main():
    # 创建输出目录
    out_dir1 = 'frame1'
    out_dir2 = 'frame2'
    out_dir3 = 'frame3'
    os.makedirs(out_dir1, exist_ok=True)
    os.makedirs(out_dir2, exist_ok=True)
    os.makedirs(out_dir3, exist_ok=True)

    # 获取当前目录下所有以 .fasta 或 .fa 结尾的文件
    files = os.listdir('.')
    fasta_files = [f for f in files if f.lower().endswith('.fasta') or f.lower().endswith('.fa')]

    if not fasta_files:
        print("当前文件夹中未发现fasta文件。")
        return

    # 对每个fasta文件进行处理
    for fasta_file in fasta_files:
        process_fasta_file(fasta_file, out_dir1, out_dir2, out_dir3)
        print(f"已处理文件：{fasta_file}")

if __name__ == "__main__":
    main()