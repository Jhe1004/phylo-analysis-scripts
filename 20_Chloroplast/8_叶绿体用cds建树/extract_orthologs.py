# -*- coding: utf-8 -*-
"""
extract_orthologs.py (Python 3.6+ Compatible, Final Robust Version)

该脚本用于从多个GenBank格式的植物叶绿体基因组文件中提取直系同源基因。

最终修复与功能:
- 智能判断基因类型，优先级: CDS > rRNA > tRNA > gene。
- [最终修复] 采用全新的、更稳健的序列提取逻辑：
  - 1. 使用 BLAST 进行初步定位，找到同源基因的大致区域。
  - 2. 在目标 GenBank 注释中精确查找所有与BLAST命中区域重叠的、同名的基因特征。
  - 3. 从候选特征中，严格排除所有假基因(pseudogene)。
  - 4. 按照特征类型优先级和序列长度(最长优先)进行双重排序，选择唯一最佳特征进行提取。
  - 此方法可完美处理含内含子(join)、反向链(complement)、以及ycf1/ndh等具有多个或模糊注释的复杂情况，贯彻“宁缺毋滥”原则。
- 文件名清理函数保留小数点（如 rrn4.5）。
- 自动选择参考基因组（如果未指定），兼容 .gb 和 .gbk 文件。
"""
import os
import sys
import argparse
import subprocess
import shutil
import re
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML
except ImportError:
    print("错误：Biopython 库未安装。请使用 'pip install biopython' 命令进行安装。")
    sys.exit(1)


def check_blast_installed():
    """检查 NCBI BLAST+ 是否已安装并可用。"""
    try:
        subprocess.run(['makeblastdb', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True)
        subprocess.run(['blastn', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("错误：找不到 NCBI BLAST+ 命令 ('makeblastdb' 或 'blastn')。")
        print("请确保您已经成功安装了 NCBI BLAST+，并将其添加到了系统的 PATH 环境变量中。")
        return False

def sanitize_filename(name):
    """清理字符串，使其成为一个合法的基因文件名（保留点）。"""
    return "".join(c for c in name if c.isalnum() or c in ('_', '-', '.')).rstrip()

def sanitize_sample_name(name):
    """清理样本名，将不安全的字符替换为下划线。"""
    name = re.sub(r'[^a-zA-Z0-9_-]', '_', name)
    name = re.sub(r'__+', '_', name)
    name = name.strip('_')
    return name

def create_blast_db(gb_files, temp_dir):
    """为所有 GenBank 文件创建 BLAST 数据库。"""
    db_paths = {}
    print("\n--- 步骤 1: 为所有 GenBank 文件创建 BLAST 数据库 ---")
    for gb_file in gb_files:
        original_stem = Path(gb_file).stem
        sample_name = sanitize_sample_name(original_stem)
        if original_stem != sample_name:
            print(f"  - 文件名清理: '{original_stem}' -> '{sample_name}'")
        fasta_path = temp_dir / f"{sample_name}.fasta"
        db_path = temp_dir / sample_name
        try:
            SeqIO.convert(gb_file, 'genbank', fasta_path, 'fasta')
            cmd = ['makeblastdb', '-in', str(fasta_path), '-dbtype', 'nucl', '-out', str(db_path)]
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            db_paths[sample_name] = str(db_path)
            print(f"  [成功] 为 {sample_name} 创建了数据库。")
        except Exception as e:
            print(f"  [失败] 为 {sample_name} 创建数据库时出错: {e}")
    return db_paths

def get_best_feature(record, gene_name, candidates=None):
    """
    [最终版逻辑] 从记录或候选列表中为给定基因名选择最佳特征。
    1. 过滤掉所有伪基因。
    2. 按类型优先级 (CDS > rRNA > tRNA > gene) 和长度 (最长优先) 排序。
    3. 返回最优选。
    """
    if candidates is None:
        candidates = [f for f in record.features if 'gene' in f.qualifiers and f.qualifiers['gene'][0] == gene_name]
    
    # 1. 严格排除伪基因
    valid_candidates = [
        f for f in candidates 
        if 'pseudo' not in f.qualifiers and 'pseudogene' not in f.qualifiers
    ]

    if not valid_candidates:
        return None, None # 宁缺毋滥

    # 2. 按优先级和长度排序
    priority_order = ['CDS', 'rRNA', 'tRNA', 'gene']
    
    def sort_key(feature):
        try:
            type_priority = priority_order.index(feature.type)
        except ValueError:
            type_priority = len(priority_order) # 未知类型优先级最低
        
        # 长度作为第二排序标准，越长越优先
        length = len(feature)
        return (type_priority, -length) # 负号表示降序

    valid_candidates.sort(key=sort_key)
    
    best_feature = valid_candidates[0]
    extraction_type = best_feature.type
    
    return best_feature, extraction_type

def find_ortholog_feature_in_blast_hit(target_record, gene_name_raw, hit_start, hit_end):
    """在BLAST命中区域内查找并提取最合适的基因特征。"""
    hit_location = set(range(hit_start, hit_end + 1))
    
    overlapping_candidates = []
    for feature in target_record.features:
        if 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] == gene_name_raw:
            feature_location = set(range(int(feature.location.start), int(feature.location.end)))
            if hit_location.intersection(feature_location):
                overlapping_candidates.append(feature)

    if not overlapping_candidates:
        return None
        
    best_feature, _ = get_best_feature(target_record, gene_name_raw, candidates=overlapping_candidates)
    
    if best_feature:
        return best_feature.extract(target_record.seq)
    else:
        return None

def main(args):
    """主执行函数"""
    if not check_blast_installed():
        sys.exit(1)

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    evalue_threshold = float(args.evalue)

    if not input_dir.is_dir():
        print(f"错误：输入目录 '{input_dir}' 不存在。")
        sys.exit(1)

    gb_files = sorted(list(input_dir.glob('*.gb'))) + sorted(list(input_dir.glob('*.gbk')))
    if not gb_files:
        print(f"错误：在目录 '{input_dir}' 中没有找到 .gb 或 .gbk 文件。")
        sys.exit(1)

    reference_file_path = None
    if args.reference:
        reference_file_path = Path(args.reference)
        if not reference_file_path.is_file():
            reference_file_path = input_dir / reference_file_path
            if not reference_file_path.is_file():
                print(f"错误：指定的参考文件 '{args.reference}' 不存在。")
                sys.exit(1)
    else:
        reference_file_path = gb_files[0]
        print(f"\n提示：未指定参考基因组，已自动选择 '{reference_file_path.name}' 作为参考。")

    output_dir.mkdir(exist_ok=True)
    temp_dir = output_dir / "temp_blast_db"
    if temp_dir.exists(): shutil.rmtree(temp_dir)
    temp_dir.mkdir(exist_ok=True)

    db_paths = create_blast_db(gb_files, temp_dir)
    if not db_paths:
        shutil.rmtree(temp_dir)
        sys.exit("错误：未能成功创建任何 BLAST 数据库，脚本终止。")

    print("\n--- 步骤 2: 读取所有 GenBank 文件记录 ---")
    all_records = {}
    for gb_file in gb_files:
        original_stem = Path(gb_file).stem
        sample_name = sanitize_sample_name(original_stem)
        try:
            record = SeqIO.read(gb_file, "genbank")
            all_records[sample_name] = record
            print(f"  [成功] 已读取 {sample_name} 的序列记录。")
        except Exception as e:
            print(f"  [失败] 读取文件 {gb_file} 时出错: {e}")
            if sample_name in db_paths:
                del db_paths[sample_name]

    print("\n--- 步骤 3: 开始提取直系同源基因 ---")
    original_ref_stem = Path(reference_file_path).stem
    reference_name = sanitize_sample_name(original_ref_stem)
    reference_record = all_records.get(reference_name)

    if not reference_record:
        shutil.rmtree(temp_dir)
        sys.exit(f"错误：无法从 all_records 字典中找到参考记录 '{reference_name}'。")
    
    all_gene_names = sorted(list(set(f.qualifiers['gene'][0] for f in reference_record.features if 'gene' in f.qualifiers)))
    
    print(f"在参考文件 '{reference_file_path.name}' 中找到 {len(all_gene_names)} 个唯一的基因进行提取。")

    for gene_name_raw in all_gene_names:
        
        feature_to_extract, extraction_type = get_best_feature(reference_record, gene_name_raw)

        if not feature_to_extract:
            print(f"\n正在处理基因: {gene_name_raw}")
            print(f"  - 警告: 在参考基因组中未找到 '{gene_name_raw}' 的有效特征(非伪基因)，已跳过。")
            continue
        
        gene_name = sanitize_filename(gene_name_raw)
        if not gene_name:
            print(f"  - 警告: 基因名 '{gene_name_raw}' 清理后为空，已跳过。")
            continue
        
        print(f"\n正在处理基因: {gene_name_raw} (提取类型: {extraction_type})")

        ref_gene_seq = feature_to_extract.extract(reference_record.seq)
        
        query_fasta = temp_dir / f"query_{gene_name}.fasta"
        with open(query_fasta, "w") as f:
            f.write(f">{gene_name}_ref\n{ref_gene_seq}\n")

        orthologs = [SeqRecord(ref_gene_seq, id=reference_name, description=f"gene={gene_name_raw} type={extraction_type}")]

        for sample_name, db_path in db_paths.items():
            if sample_name == reference_name:
                continue

            target_record = all_records.get(sample_name)
            if not target_record:
                print(f"  - 警告: 找不到样本 {sample_name} 的序列记录，跳过。")
                continue

            blast_result_xml = temp_dir / f"result_{sample_name}.xml"
            blastn_cline = NcbiblastnCommandline(query=str(query_fasta), db=db_path,
                                                 evalue=evalue_threshold,
                                                 outfmt=5, out=str(blast_result_xml))
            
            stdout, stderr = blastn_cline()

            try:
                best_hit_seq = None
                with open(blast_result_xml) as result_handle:
                    blast_records = NCBIXML.parse(result_handle)
                    blast_record = next(blast_records, None)

                    if blast_record and blast_record.alignments:
                        hsp = blast_record.alignments[0].hsps[0]

                        if hsp.expect <= evalue_threshold:
                            start, end = min(hsp.sbjct_start, hsp.sbjct_end), max(hsp.sbjct_start, hsp.sbjct_end)
                            extracted_seq = find_ortholog_feature_in_blast_hit(target_record, gene_name_raw, start, end)
                            
                            if extracted_seq:
                                best_hit_seq = extracted_seq
                            else:
                                print(f"  - 在 {sample_name} 中未找到对应的有效基因特征。")
                        else:
                             print(f"  - {sample_name}: 最佳匹配的 E-value ({hsp.expect:.2e}) 高于阈值，跳过。")
                
                if best_hit_seq:
                    ortholog_record = SeqRecord(best_hit_seq, id=sample_name, description=f"ortholog_of={gene_name_raw}")
                    orthologs.append(ortholog_record)
                    print(f"  - 在 {sample_name} 中找到并提取了完整的同源基因。")
                # else: # 不再打印 "未找到满足条件的" 警告，因为上面已经打印了更具体的信息
                #     pass

            except (ValueError, StopIteration):
                print(f"  - 在 {sample_name} 中未找到任何 BLAST 匹配项。")
            except Exception as e:
                print(f"  - 错误: 解析 {sample_name} 的 BLAST 结果时出错: {e}")

        if len(orthologs) > 1:
            output_fasta_path = output_dir / f"{gene_name}.fasta"
            SeqIO.write(orthologs, output_fasta_path, "fasta")
            print(f"  [完成] {len(orthologs)} 个同源序列已保存到 {output_fasta_path}")
        else:
            print(f"  [跳过] 未能为基因 {gene_name_raw} 找到任何其他同源序列。")

    print("\n--- 步骤 4: 清理临时文件和数据库 ---")
    try:
        shutil.rmtree(temp_dir)
        print("  临时文件已成功删除。")
    except OSError as e:
        print(f"  错误: 删除临时目录 {temp_dir} 失败: {e}")

    print("\n脚本执行完毕！")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="从多个 GenBank 文件中提取直系同源基因 (Python 3.6+ 兼容, 优先提取CDS)。",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
使用示例:
  # 自动选择参考文件
  python extract_orthologs.py -i ./input_data -o ./output

  # 手动指定参考文件
  python extract_orthologs.py -i ./input_data -o ./output -r my_favorite_genome.gb
"""
    )
    parser.add_argument('-r', '--reference', type=str, required=False, 
                        help='参考基因组的 GenBank 文件名。如果省略，将自动选择输入目录中的第一个文件作为参考。')
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='包含所有 GenBank 文件的输入目录。')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='存放输出 FASTA 文件的目录。')
    parser.add_argument('-e', '--evalue', type=str, default='1e-10', help='BLASTn 搜索的 E-value 阈值 (默认: 1e-10)。')
    
    args = parser.parse_args()
    main(args)