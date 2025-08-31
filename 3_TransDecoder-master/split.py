import os
from Bio import SeqIO

def split_sequence(sequence, max_length=40000):
    """将序列分割为指定长度的子序列."""
    return [sequence[i:i + max_length] for i in range(0, len(sequence), max_length)]

def split_fasta(input_fasta, output_fasta, max_length=40000):
    """读取输入的FASTA文件并将长序列分割后写入新的FASTA文件."""
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if len(record.seq) > max_length:
                sub_sequences = split_sequence(record.seq, max_length)
                for i, sub_seq in enumerate(sub_sequences):
                    new_record = record[:]
                    new_record.id = f"{record.id}_part{i+1}"
                    new_record.seq = sub_seq
                    SeqIO.write(new_record, output_handle, "fasta")
            else:
                SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    # 获取当前目录下所有.fasta文件
    fasta_files = [f for f in os.listdir('.') if f.endswith('.fasta')]

    for fasta_file in fasta_files:
        output_fasta = fasta_file.replace(".fasta", ".fas")  # 为输出文件命名
        split_fasta(fasta_file, output_fasta)
        print(f"已处理文件: {fasta_file}, 结果保存在: {output_fasta}")