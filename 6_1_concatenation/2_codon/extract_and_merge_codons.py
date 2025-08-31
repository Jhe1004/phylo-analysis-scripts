from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

if len(sys.argv) < 2:
    print("用法: python extract_and_merge_codons.py input_alignment.fasta")
    sys.exit(1)

input_file = sys.argv[1]

# 读取输入对齐文件
records = list(SeqIO.parse(input_file, "fasta"))
if not records:
    print("输入文件中没有序列记录，结束。")
    sys.exit(1)

# 分别存储1,2,3位密码子位点序列
codon1_records = []
codon2_records = []
codon3_records = []

for rec in records:
    seq_str = str(rec.seq)
    codon1_seq = ''.join(seq_str[i] for i in range(0, len(seq_str), 3) if i < len(seq_str))
    codon2_seq = ''.join(seq_str[i] for i in range(1, len(seq_str), 3) if i < len(seq_str))
    codon3_seq = ''.join(seq_str[i] for i in range(2, len(seq_str), 3) if i < len(seq_str))

    codon1_records.append(SeqRecord(Seq(codon1_seq), id=rec.id, description="codon1_positions"))
    codon2_records.append(SeqRecord(Seq(codon2_seq), id=rec.id, description="codon2_positions"))
    codon3_records.append(SeqRecord(Seq(codon3_seq), id=rec.id, description="codon3_positions"))

# 输出三个新的对齐文件
SeqIO.write(codon1_records, "codon1_aln.fasta", "fasta")
SeqIO.write(codon2_records, "codon2_aln.fasta", "fasta")
SeqIO.write(codon3_records, "codon3_aln.fasta", "fasta")

print("已生成 codon1_aln.fasta, codon2_aln.fasta, codon3_aln.fasta 三个文件。")

# 现在将三个对齐文件合并为一个
# 合并规则：对于同一物种(ID)的序列，将codon1, codon2, codon3三段接在一起

# 假设三个文件的物种ID和顺序是一致的（因为来自同一个源对齐）
# 在实际情况中，建议确认三个文件的序列ID顺序一致性
codon1_recs = list(SeqIO.parse("codon1_aln.fasta", "fasta"))
codon2_recs = list(SeqIO.parse("codon2_aln.fasta", "fasta"))
codon3_recs = list(SeqIO.parse("codon3_aln.fasta", "fasta"))

if len(codon1_recs) != len(codon2_recs) or len(codon2_recs) != len(codon3_recs):
    print("警告：三个codon对齐文件的序列数不匹配，请检查数据一致性。")
    sys.exit(1)

# 将三段序列合并
combined_records = []
for r1, r2, r3 in zip(codon1_recs, codon2_recs, codon3_recs):
    combined_seq = str(r1.seq) + str(r2.seq) + str(r3.seq)
    combined_records.append(SeqRecord(Seq(combined_seq), id=r1.id, description="combined_codon"))

# 写出合并后的文件
SeqIO.write(combined_records, "combined_aln.fasta", "fasta")
print("已生成合并的对齐文件：combined_aln.fasta")

# 生成partitions.txt文件
# 计算各个区段的长度
if len(combined_records) > 0:
    total_length = len(combined_records[0].seq)
    codon1_length = len(codon1_recs[0].seq)
    codon2_length = len(codon2_recs[0].seq)
    codon3_length = len(codon3_recs[0].seq)

    # 写入partitions.txt
    with open("partitions.txt", "w") as f:
        f.write(f"DNA, codon1 = 1-{codon1_length}\n")
        f.write(f"DNA, codon2 = {codon1_length+1}-{codon1_length+codon2_length}\n")
        f.write(f"DNA, codon3 = {codon1_length+codon2_length+1}-{codon1_length+codon2_length+codon3_length}\n")

    print("已生成 RAxML 所需的 partitions.txt 文件。")

else:
    print("合并文件中无序列，未生成partitions.txt")