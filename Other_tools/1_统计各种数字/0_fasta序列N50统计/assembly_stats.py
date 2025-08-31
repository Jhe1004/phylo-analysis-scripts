import os
from Bio import SeqIO

def calculate_n50(lengths):
    """计算N50"""
    lengths_sorted = sorted(lengths, reverse=True)
    total_length = sum(lengths)
    half_total_length = total_length / 2
    running_sum = 0
    for length in lengths_sorted:
        running_sum += length
        if running_sum >= half_total_length:
            return length
    return None

def process_fasta_file(fasta_path):
    """统计FASTA文件中的序列数目和N50"""
    seq_lengths = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_lengths.append(len(record.seq))
    sequence_count = len(seq_lengths)
    n50 = calculate_n50(seq_lengths)
    return sequence_count, n50

def scan_directory_for_fasta(directory):
    """扫描目录并统计每个FASTA文件中的序列数和N50"""
    results = []
    for file_name in os.listdir(directory):
        if file_name.endswith(".fasta") or file_name.endswith(".fa") or file_name.endswith(".cds"):
            fasta_path = os.path.join(directory, file_name)
            seq_count, n50 = process_fasta_file(fasta_path)
            results.append((file_name, seq_count, n50))
    return results

def main(directory):
    results = scan_directory_for_fasta(directory)
    if results:
        print(f"{'File':<30} {'Sequences':<10} {'N50':<10}")
        print("-" * 50)
        for file_name, seq_count, n50 in results:
            print(f"{file_name:<30} {seq_count:<10} {n50:<10}")
    else:
        print("没有找到FASTA文件")

# 使用你的当前目录
if __name__ == "__main__":
    fasta_directory = '.'  # 当前目录，或指定其他目录
    main(fasta_directory)