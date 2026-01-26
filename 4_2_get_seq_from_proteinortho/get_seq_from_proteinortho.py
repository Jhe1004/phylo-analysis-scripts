# -*- coding: utf-8 -*-
import os          
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

# 内置参数
proteinortho_output_file = "myproject.proteinortho.tsv"  # 输入的proteinortho文件
output_seq_dir = "output2"  # 输出文件夹
output_summary_file = "selected_orthologs.tsv" # 新增参数：输出的筛选结果文件名
max_allow_missing_species_proportion = 0.000001  # 最大缺失物种比例
max_allow_low_copy_num = 1  # 最大低拷贝数
min_alg_conn = 0.5  # 新增参数：Alg.-Conn.的最小值，范围0-1

# make output dir
try:
    os.makedirs(os.getcwd() + os.sep + output_seq_dir)
except FileExistsError:
    print(f"输出文件夹 \"{output_seq_dir}\" 已存在，脚本将继续并在其中写入文件。")
except Exception as e:
    print(f"创建输出文件夹时出错: {e}")
    sys.exit()


# Import the FASTA sequence into memory
print("正在将序列文件加载到内存中...")
file_dict = {}  
for each in os.listdir(os.getcwd()):        
    if ".pep" in each:
        try:
            file_dict[each] = SeqIO.to_dict(SeqIO.parse(each, "fasta"))
        except Exception as e:
            print(f"无法解析文件 {each}: {e}")
    elif ".cds" in each:
        try:
            file_dict[each] = SeqIO.to_dict(SeqIO.parse(each, "fasta"))
        except Exception as e:
            print(f"无法解析文件 {each}: {e}")

if not file_dict:
    print("错误：在当前目录下未找到任何 .pep 或 .cds 文件。")
    sys.exit()
print("序列文件加载完毕。")

# get sequences from ".cds" or ".pep" files.            
def get_nucleotide(IDs): 
    try:
        cds_file = IDs.split("++")[0][:-3] + "cds"
        seq_id = IDs.split("++")[1]
        return str(file_dict[cds_file][seq_id].seq)
    except KeyError:
        # print(f"警告：在文件 {cds_file} 中找不到核苷酸序列ID {seq_id}")
        return ""

def get_protein(IDs): 
    try:
        pep_file = IDs.split("++")[0][:-3] + "pep"
        seq_id = IDs.split("++")[1]
        return str(file_dict[pep_file][seq_id].seq)
    except KeyError:
        # print(f"警告：在文件 {pep_file} 中找不到蛋白质序列ID {seq_id}")
        return ""

# main analysis
try:
    df = pd.read_csv(proteinortho_output_file, sep='\t', header=0)
except FileNotFoundError:
    print(f"错误：找不到输入的proteinortho文件 '{proteinortho_output_file}'")
    sys.exit()

columns_values = df.columns.values
total_num = len(columns_values) - 3

print("开始处理proteinortho文件并提取单拷贝直系同源基因...")
processed_count = 0
selected_rows = [] # 用于存储所有符合条件的行

for i in range(len(df)):
    row = df.iloc[i,:]
    
    # 筛选条件1：缺失物种比例
    if (row[3:].tolist().count("*"))/len(row[3:]) <= max_allow_missing_species_proportion:
        
        # 筛选条件2：Alg.-Conn. 值
        if row['Alg.-Conn.'] >= min_alg_conn:
            
            nucleotide_fasta_path = os.path.join(os.getcwd(), output_seq_dir, f"ortho{i}_cds.fasta")
            protein_fasta_path = os.path.join(os.getcwd(), output_seq_dir, f"ortho{i}_pep.fasta")

            with open(nucleotide_fasta_path, "w") as nucleotide_fasta, open(protein_fasta_path, "w") as protein_fasta:
                valid_ortholog = True
                sequences_to_write = []

                for n in range(3, len(row)):
                    if row.iloc[n] != "*":
                        if len(row.iloc[n].split(",")) == 1:
                            sequences_to_write.append( (columns_values[n], row.iloc[n]) )

                        elif len(row.iloc[n].split(",")) > max_allow_low_copy_num:
                            valid_ortholog = False
                            break
                        else:
                            a = ""
                            b = ""
                            for each_seq in row.iloc[n].split(","):
                                seq = get_nucleotide(columns_values[n] + "++" + each_seq)
                                if len(seq) > len(a):
                                    a = seq
                                    b = each_seq
                            sequences_to_write.append( (columns_values[n], b) )
                
                if valid_ortholog:
                    processed_count += 1
                    selected_rows.append(row) # 将符合条件的行添加到列表中
                    for species_name, seq_id in sequences_to_write:
                        full_id = f"{species_name}++{seq_id}"
                        nucleotide_seq = get_nucleotide(full_id)
                        protein_seq = get_protein(full_id)
                        
                        if nucleotide_seq:
                            nucleotide_fasta.write(f">{species_name}\n")
                            nucleotide_fasta.write(f"{nucleotide_seq}\n")
                        
                        if protein_seq:
                            protein_fasta.write(f">{species_name}\n")
                            protein_fasta.write(f"{protein_seq}\n")
                else:
                    # 如果因为拷贝数过多而无效，则删除已创建的文件
                    pass # with open 会自动关闭文件，之后可以安全删除

            if not valid_ortholog:
                # 确保在验证失败时删除空文件
                if os.path.exists(nucleotide_fasta_path):
                    os.remove(nucleotide_fasta_path)
                if os.path.exists(protein_fasta_path):
                    os.remove(protein_fasta_path)

# 在循环结束后，将筛选出的行写入新的tsv文件
if selected_rows:
    selected_df = pd.DataFrame(selected_rows)
    summary_file_path = os.path.join(os.getcwd(), output_seq_dir, output_summary_file)
    selected_df.to_csv(summary_file_path, sep='\t', index=False)
    print(f"已将筛选出的 {len(selected_rows)} 行记录保存到文件: {summary_file_path}")
else:
    print("没有找到符合条件的直系同源基因簇。")


print(f"处理完成！共提取了 {processed_count} 个符合条件的直系同源基因簇。")

