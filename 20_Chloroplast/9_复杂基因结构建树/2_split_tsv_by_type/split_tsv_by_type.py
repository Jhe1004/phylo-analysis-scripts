#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path


INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/9_复杂基因结构建树/2_split_tsv_by_type/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/9_复杂基因结构建树/2_split_tsv_by_type/output"
INPUT_TSV = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/9_复杂基因结构建树/2_split_tsv_by_type/input/myproject.proteinortho.tsv"
OUTPUT_CDS = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/9_复杂基因结构建树/2_split_tsv_by_type/output/proteinortho_cds.tsv"
OUTPUT_RNA = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/9_复杂基因结构建树/2_split_tsv_by_type/output/proteinortho_rna.tsv"
OUTPUT_IGS = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/20_Chloroplast/9_复杂基因结构建树/2_split_tsv_by_type/output/proteinortho_igs.tsv"

def get_fragment_type(fragment_id):
    """
    根据序列ID中的关键词推断其类型。
    rRNA 和 tRNA 将被合并为 'rna'。
    """
    # 检查的顺序很重要
    
    # 1. 检查 IGS (最明确)
    if "_IGS_" in fragment_id:
        return "igs"
    
    # 2. 检查 rRNA (例如: _rrn16S_)
    # 使用 .lower() 确保 _rrn, _RRN 都匹配
    if "_rrn" in fragment_id.lower():
        return "rna"  # <-- 归类为 "rna"
        
    # 3. 检查 tRNA (例如: _trnK-UUU_)
    if "_trn" in fragment_id.lower():
        return "rna"  # <-- 归类为 "rna"
        
    # 4. 如果都不是，则假定为 CDS
    # (例如: _atpE_, _rps4_, _rpoC1_exon1_)
    return "cds"

# --- 主函数 ---

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_PATH = SCRIPT_DIR / INPUT_DIRECTORY / INPUT_TSV
OUTPUT_DIR = SCRIPT_DIR / OUTPUT_DIRECTORY


def split_tsv():
    """
    主函数：读取TSV文件，将其拆分为 CDS, RNA (rRNA+tRNA), 和 IGS 三个文件。
    """
    print(f"开始处理: {INPUT_PATH}...")
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"未找到输入文件: {INPUT_PATH}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        
    # 存储标题行
    header_lines = []
    
    # 为每种类型创建一个列表来存储对应的行
    lines_data = {
        "cds": [],
        "rna": [],  # <-- 合并后的类别
        "igs": [],
        "unknown": [] # 存储所有物种均为 '*' 的行
    }
    
    try:
        # --- 步骤 1: 读取和分类所有行 ---
        with INPUT_PATH.open('r', encoding='utf-8') as f_in:
            for line in f_in:
                # 1. 存储注释/标题行
                if line.startswith('#'):
                    header_lines.append(line)
                    continue 
                
                # 跳过空行
                if not line.strip():
                    continue

                # 2. 处理数据行
                cols = line.strip().split('\t')
                group_type = "unknown"
                
                # 从第4列 (索引3) 开始查找第一个非'*'的单元格
                for cell in cols[3:]:
                    if cell != '*':
                        # 找到了第一个有效的片段ID
                        first_fragment_id = cell.split(',')[0]
                        group_type = get_fragment_type(first_fragment_id)
                        break # 已经确定了该行(OG)的类型，跳出循环
                
                # 3. 将整行添加到对应的列表中
                lines_data[group_type].append(line)

        # --- 步骤 2: 将分类后的行写入各自的文件 ---
        
        # 写入 CDS 文件
        with (OUTPUT_DIR / OUTPUT_CDS).open('w', encoding='utf-8') as f_out:
            f_out.writelines(header_lines) # 写入表头
            f_out.writelines(lines_data["cds"]) # 写入数据
            
        # 写入 RNA 文件
        with (OUTPUT_DIR / OUTPUT_RNA).open('w', encoding='utf-8') as f_out:
            f_out.writelines(header_lines)
            f_out.writelines(lines_data["rna"])

        # 写入 IGS 文件
        with (OUTPUT_DIR / OUTPUT_IGS).open('w', encoding='utf-8') as f_out:
            f_out.writelines(header_lines)
            f_out.writelines(lines_data["igs"])

        print("\n处理完成。文件已成功拆分：")
        print(f" - {OUTPUT_CDS} ({len(lines_data['cds'])} 个同源组)")
        print(f" - {OUTPUT_RNA} ({len(lines_data['rna'])} 个同源组)")
        print(f" - {OUTPUT_IGS} ({len(lines_data['igs'])} 个同源组)")

        if lines_data["unknown"]:
            print(f"\n警告: 有 {len(lines_data['unknown'])} 行无法确定类型 (可能所有物种均为'*')。")

    except Exception as e:
        print(f"发生错误: {e}")

# --- 运行脚本 ---
if __name__ == "__main__":
    split_tsv()
