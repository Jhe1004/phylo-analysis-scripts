# -*- coding: utf-8 -*-
"""
extract_orthologs.py (Python 3.6+ Compatible, with CDS extraction, XML fix, Filename Sanitize, Auto-Reference)

该脚本用于从多个GenBank格式的植物叶绿体基因组文件中提取直系同源基因。

新功能：
- 脚本会智能判断基因类型。
- 如果是蛋白质编码基因，则优先提取其编码序列 (CDS, 即拼接后的外显子)。
- 如果是非编码基因 (如 tRNA, rRNA) 或没有CDS注释，则提取基因全长。
- 修复了因提前中断XML解析而产生的ExpatError警告。
- 新增：自动清理样本文件名，将空格、点等特殊字符替换为下划线，以增加健壮性。
- 新增：如果未指定参考基因组，则自动选择输入目录中的第一个 GenBank 文件作为参考。

工作原理：
1. 指定一个样本作为参考基因组（若不指定则自动选择）。
2. 脚本为当前文件夹中所有的 .gb 文件创建本地 BLAST 数据库。
3. 遍历参考基因组中的每一个基因，并判断是否提取CDS或基因全长。
4. 将该序列作为 query，使用 blastn 在其他所有样本的数据库中搜索同源序列。
5. 提取每个样本中比对结果最好 (best hit) 的序列。
6. 将所有样本中找到的该基因序列汇总，保存到一个以该基因命名的 FASTA 文件中。
7. 对参考基因组中的所有基因重复此过程。
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
        subprocess.run(
            ['makeblastdb', '-version'], 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
            check=True, universal_newlines=True
        )
        subprocess.run(
            ['blastn', '-version'], 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
            check=True, universal_newlines=True
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("错误：找不到 NCBI BLAST+ 命令 ('makeblastdb' 或 'blastn')。")
        print("请确保您已经成功安装了 NCBI BLAST+，并将其添加到了系统的 PATH 环境变量中。")
        print("对于 Conda 用户，推荐使用 'conda install -c bioconda blast' 命令安装。")
        return False

def sanitize_filename(name):
    """清理字符串，使其成为一个合法的基因文件名。"""
    return "".join(c for c in name if c.isalnum() or c in ('_', '-')).rstrip()

def sanitize_sample_name(name):
    """
    清理样本名，将不安全的字符替换为下划线。
    - 替换所有非字母数字、非下划线、非连字符的字符为 '_'
    - 将多个连续的 '_' 合并为一个
    - 移除开头和结尾的 '_'
    """
    name = re.sub(r'[^a-zA-Z0-9_-]', '_', name)
    name = re.sub(r'__+', '_', name)
    name = name.strip('_')
    return name

def create_blast_db(gb_files, temp_dir):
    """为所有 GenBank 文件创建 BLAST 数据库。"""
    db_paths = {}
    print("\n--- 步骤 1: 为所有 GenBank 文件创建 BLAST 数据库 ---")
    print("注意：文件名中的特殊字符 (如空格, 点) 将被替换为下划线 '_'。")
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
        except subprocess.CalledProcessError as e:
            print(f"  [失败] 为 {sample_name} 创建数据库时出错。")
            print(f"  命令: {' '.join(e.cmd)}")
            print(f"  错误输出 (stderr):\n{e.stderr}")
            print("  请检查 GenBank 文件格式是否正确或 BLAST+ 是否安装正确。")
        except Exception as e:
            print(f"  [失败] 为 {sample_name} 创建数据库时出现意外错误: {e}")

    return db_paths

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

    # --- 修改: 自动选择参考文件的逻辑 ---
    reference_file_path = None
    if args.reference:
        # 如果用户指定了参考文件
        reference_file_path = Path(args.reference)
        if not reference_file_path.is_file():
            reference_file_path = input_dir / reference_file_path
            if not reference_file_path.is_file():
                print(f"错误：指定的参考文件 '{args.reference}' 不存在。")
                sys.exit(1)
    else:
        # 如果用户未指定，自动选择第一个文件
        reference_file_path = gb_files[0]
        print(f"\n提示：未指定参考基因组，已自动选择 '{reference_file_path.name}' 作为参考。")
    # --- 修改结束 ---

    output_dir.mkdir(exist_ok=True)
    temp_dir = input_dir / "temp_blast_db"
    temp_dir.mkdir(exist_ok=True)

    db_paths = create_blast_db(gb_files, temp_dir)
    if not db_paths:
        print("错误：未能成功创建任何 BLAST 数据库，脚本终止。")
        shutil.rmtree(temp_dir)
        sys.exit(1)

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

    print("\n--- 步骤 3: 开始提取直系同源基因 (优先提取CDS) ---")
    original_ref_stem = Path(reference_file_path).stem
    reference_name = sanitize_sample_name(original_ref_stem)
    reference_record = all_records.get(reference_name)

    if not reference_record:
        print(f"错误：无法从 all_records 字典中找到参考记录 '{reference_name}'。")
        print(f"原始文件名: '{original_ref_stem}'")
        shutil.rmtree(temp_dir)
        sys.exit(1)
    
    gene_map = {}
    cds_map = {}
    for feature in reference_record.features:
        if 'gene' in feature.qualifiers:
            gene_name = feature.qualifiers['gene'][0]
            if feature.type == 'gene':
                gene_map[gene_name] = feature
            elif feature.type == 'CDS':
                cds_map[gene_name] = feature

    print(f"在参考文件 '{reference_file_path.name}' (使用名: {reference_name}) 中找到 {len(gene_map)} 个基因。")

    for gene_name_raw, gene_feature in sorted(gene_map.items()):
        
        feature_to_extract = None
        extraction_type = ""

        if gene_name_raw in cds_map:
            feature_to_extract = cds_map[gene_name_raw]
            extraction_type = "CDS"
        else:
            feature_to_extract = gene_feature
            extraction_type = "基因全长"
        
        gene_name = sanitize_filename(gene_name_raw)
        if not gene_name:
            print(f"  - 警告: 基因名 '{gene_name_raw}' 清理后为空，已跳过。")
            continue
        
        print(f"\n正在处理基因: {gene_name_raw} (提取 {extraction_type})")

        ref_gene_seq = feature_to_extract.extract(reference_record.seq)
        
        query_fasta = temp_dir / f"query_{gene_name}.fasta"
        with open(query_fasta, "w") as f:
            f.write(f">{gene_name}_ref\n{ref_gene_seq}\n")

        orthologs = []
        ref_seq_record = SeqRecord(ref_gene_seq, id=reference_name, description=f"gene={gene_name_raw} type={extraction_type}")
        orthologs.append(ref_seq_record)

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
                    blast_record = NCBIXML.read(result_handle)

                    if blast_record.alignments:
                        alignment = blast_record.alignments[0]
                        hsp = alignment.hsps[0]

                        if hsp.expect <= evalue_threshold:
                            start, end = hsp.sbjct_start, hsp.sbjct_end
                            if start > end:
                                start, end = end, start
                            
                            hit_seq = target_record.seq[start-1:end]

                            if hsp.strand[1] == -1:
                                hit_seq = hit_seq.reverse_complement()
                            
                            best_hit_seq = hit_seq
                        else:
                             print(f"  - {sample_name}: 最佳匹配的 E-value ({hsp.expect:.2e}) 高于阈值，已跳过。")
                
                if best_hit_seq:
                    ortholog_record = SeqRecord(best_hit_seq, id=sample_name, description=f"ortholog_of={gene_name_raw}")
                    orthologs.append(ortholog_record)
                    print(f"  - 在 {sample_name} 中找到同源基因。")
                else:
                    print(f"  - 警告: 在 {sample_name} 中未找到满足条件的同源基因。")

            except ValueError:
                print(f"  - 在 {sample_name} 中未找到任何 BLAST 匹配项。")
            except Exception as e:
                print(f"  - 错误: 解析 {sample_name} 的 BLAST 结果时出错: {e}")

        if len(orthologs) > 1:
            output_fasta_path = output_dir / f"{gene_name}.fasta"
            with open(output_fasta_path, "w") as output_handle:
                SeqIO.write(orthologs, output_handle, "fasta")
            print(f"  [完成] {len(orthologs)} 个同源序列已保存到 {output_fasta_path}")
        else:
            print(f"  [跳过] 未能为基因 {gene_name_raw} 找到任何其他同源序列，不生成FASTA文件。")


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
  python extract_orthologs.py -i ./input_genomes -o ./ortholog_results

  # 手动指定参考文件
  python extract_orthologs.py -i ./input_genomes -o ./ortholog_results -r my_favorite_genome.gb
"""
    )
    # --- 修改: reference 参数不再是必须的 ---
    parser.add_argument('-r', '--reference', type=str, required=False, 
                        help='参考基因组的 GenBank 文件名。如果省略，将自动选择输入目录中的第一个文件作为参考。')
    # --- 修改结束 ---
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='包含所有 GenBank 文件的输入目录。')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='存放输出 FASTA 文件的目录。')
    parser.add_argument('-e', '--evalue', type=str, default='1e-10', help='BLASTn 搜索的 E-value 阈值 (默认: 1e-10)。')
    
    args = parser.parse_args()
    main(args)