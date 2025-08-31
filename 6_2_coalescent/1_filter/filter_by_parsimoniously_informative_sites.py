import os
import shutil
from Bio import SeqIO

# 设置简约信息位点数量的阈值
THRESHOLD = 1000

# 设置输出文件夹名称
OUTPUT_DIR = "PI1000"

# 定义有效碱基集合
valid_bases = set(["A", "T", "C", "G", "a", "t", "c", "g"])

# 确保输出目录存在，不存在则创建
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# 获取当前目录下的所有fasta文件
fasta_files = [f for f in os.listdir('.') if f.lower().endswith(('.fasta', '.fa', '.fas'))]

for fasta_file in fasta_files:
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        # 如果fasta文件中没有序列记录，跳过
        continue
    
    seq_length = len(records[0].seq)
    
    # 检查序列长度一致性
    for rec in records:
        if len(rec.seq) != seq_length:
            print(f"警告：{fasta_file} 中的序列长度不一致")
            # 此处不终止，继续执行，但可根据实际情况进行raise异常
    
    parsimonious_sites_count = 0
    
    # 统计简约信息位点数目
    for i in range(seq_length):
        base_counts = {}
        for rec in records:
            base = str(rec.seq[i]).upper()
            if base in valid_bases:
                base_counts[base] = base_counts.get(base, 0) + 1
        
        # 判断是否为简约信息位点
        if len(base_counts) >= 2:
            frequent_bases = [b for b, count in base_counts.items() if count >= 2]
            if len(frequent_bases) >= 2:
                parsimonious_sites_count += 1
    
    print(f"{fasta_file} 的简约信息位点数目: {parsimonious_sites_count}")
    
    # 若简约信息位点数目超过阈值，则复制fasta和对应的树文件至新目录
    if parsimonious_sites_count > THRESHOLD:
        # 推断树文件名称
        base_name = os.path.splitext(fasta_file)[0]  # 去掉.fas/.fasta/.fa后缀
        tree_file = "RAxML_bipartitions." + os.path.basename(base_name)
        
        # 检查对应的树文件是否存在
        if not os.path.exists(tree_file):
            print(f"警告：未找到 {fasta_file} 对应的树文件 {tree_file}")
            # 根据需要决定是否继续，此处选择继续
        
        # 复制fasta文件至指定目录
        shutil.copy(fasta_file, OUTPUT_DIR)
        # 如果树文件存在则复制
        if os.path.exists(tree_file):
            shutil.copy(tree_file, OUTPUT_DIR)