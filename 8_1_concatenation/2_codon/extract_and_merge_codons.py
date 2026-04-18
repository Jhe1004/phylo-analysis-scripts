import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_1_concatenation/2_codon/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_1_concatenation/2_codon/output"
INPUT_FILENAME = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_1_concatenation/2_codon/input/input_alignment.fasta"
DRY_RUN = False
# ============================================================


CODON1_FILENAME = "codon1_aln.fasta"
CODON2_FILENAME = "codon2_aln.fasta"
CODON3_FILENAME = "codon3_aln.fasta"
COMBINED_FILENAME = "combined_aln.fasta"
PARTITION_FILENAME = "partitions.txt"


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    input_dir = INPUT_DIRECTORY
    output_dir = OUTPUT_DIRECTORY
    input_file = INPUT_FILENAME
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.isfile(input_file):
        print(f"错误: 输入文件不存在: {input_file}")
        sys.exit(1)

    if DRY_RUN:
        print(f"[DRY RUN] Would process: {input_file}")
        sys.exit(0)

    records = list(SeqIO.parse(input_file, "fasta"))
    if not records:
        print("输入文件中没有序列记录。")
        sys.exit(1)

    codon1_records = []
    codon2_records = []
    codon3_records = []
    for rec in records:
        seq_str = str(rec.seq)
        codon1_seq = "".join(seq_str[i] for i in range(0, len(seq_str), 3))
        codon2_seq = "".join(seq_str[i] for i in range(1, len(seq_str), 3))
        codon3_seq = "".join(seq_str[i] for i in range(2, len(seq_str), 3))
        codon1_records.append(SeqRecord(Seq(codon1_seq), id=rec.id, description="codon1_positions"))
        codon2_records.append(SeqRecord(Seq(codon2_seq), id=rec.id, description="codon2_positions"))
        codon3_records.append(SeqRecord(Seq(codon3_seq), id=rec.id, description="codon3_positions"))

    codon1_path = os.path.join(output_dir, CODON1_FILENAME)
    codon2_path = os.path.join(output_dir, CODON2_FILENAME)
    codon3_path = os.path.join(output_dir, CODON3_FILENAME)
    SeqIO.write(codon1_records, codon1_path, "fasta")
    SeqIO.write(codon2_records, codon2_path, "fasta")
    SeqIO.write(codon3_records, codon3_path, "fasta")

    combined_records = []
    for r1, r2, r3 in zip(codon1_records, codon2_records, codon3_records):
        combined_seq = str(r1.seq) + str(r2.seq) + str(r3.seq)
        combined_records.append(SeqRecord(Seq(combined_seq), id=r1.id, description="combined_codon"))

    combined_path = os.path.join(output_dir, COMBINED_FILENAME)
    SeqIO.write(combined_records, combined_path, "fasta")

    codon1_length = len(codon1_records[0].seq)
    codon2_length = len(codon2_records[0].seq)
    codon3_length = len(codon3_records[0].seq)
    partition_path = os.path.join(output_dir, PARTITION_FILENAME)
    with open(partition_path, "w", encoding="utf-8") as handle:
        handle.write(f"DNA, codon1 = 1-{codon1_length}\n")
        handle.write(f"DNA, codon2 = {codon1_length+1}-{codon1_length+codon2_length}\n")
        handle.write(
            f"DNA, codon3 = {codon1_length+codon2_length+1}-{codon1_length+codon2_length+codon3_length}\n"
        )

    print(f"输出完成: {output_dir}")


if __name__ == "__main__":
    main()
