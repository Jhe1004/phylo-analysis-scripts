import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def filter_sequences(input_dir):
    for filename in os.listdir(input_dir):
        if filename.endswith(".pep"):
            process_pep_fasta(os.path.join(input_dir, filename), 200, "filtered")
        elif filename.endswith(".cds"):
            process_fasta(os.path.join(input_dir, filename), 600, "filtered")

def process_fasta(filepath, min_length, suffix):
    base_name, ext = os.path.splitext(os.path.basename(filepath))
    output_filename = os.path.join(os.path.dirname(filepath), f"{base_name}_{suffix}{ext}")
    
    with open(output_filename, 'w') as output_handle:
        for record in SeqIO.parse(filepath, "fasta"):
            if len(record.seq) > min_length:
                SeqIO.write(record, output_handle, "fasta")
    print(f"Filtered sequences saved to {output_filename}")

def process_pep_fasta(filepath, min_length, suffix):
    base_name, ext = os.path.splitext(os.path.basename(filepath))
    output_filename = os.path.join(os.path.dirname(filepath), f"{base_name}_{suffix}{ext}")
    
    with open(output_filename, 'w') as output_handle:
        for record in SeqIO.parse(filepath, "fasta"):
            sequence = str(record.seq)
            # 移除末尾的'*'符号
            if sequence.endswith('*'):
                sequence = sequence[:-1]
            if len(sequence) > min_length and len(sequence) > 0:  # 确保序列有效且不为空
                # 创建新的SeqRecord对象，避免修改原始record对象
                new_record = SeqRecord(Seq(sequence), id=record.id, description=record.description)
                SeqIO.write(new_record, output_handle, "fasta")
    print(f"Filtered and cleaned sequences saved to {output_filename}")

if __name__ == "__main__":
    input_directory = "./"  # 修改为你的文件夹路径
    filter_sequences(input_directory)