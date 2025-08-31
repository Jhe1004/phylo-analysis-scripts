# -*- coding: utf-8 -*-
"""
一个用于将 FASTA 格式的多序列比对文件转换为 VCF 格式的 Python 脚本。

该脚本会：
1. 读取一个包含多条序列的比对文件（FASTA 格式）。
2. 将其中一条序列作为参考（reference）。
3. 比较其他所有序列与参考序列，找出变异位点。
4. 生成一个 VCF 格式的输出文件，其中包含每个样本在变异位点的基因型信息。

使用方法:
    python fasta_to_vcf.py -i your_alignment.fasta -o output.vcf -r reference_sequence_id

参数说明:
    -i, --input:    输入的 FASTA 文件路径。
    -o, --output:   输出的 VCF 文件路径。
    -r, --ref:      在 FASTA 文件中用作参考序列的序列 ID (例如 >sequence1 中的 "sequence1")。

依赖:
    - Biopython (请先通过 'pip install biopython' 安装)
"""

import argparse
import sys
from datetime import datetime
try:
    from Bio import SeqIO
except ImportError:
    print("错误：本脚本需要 Biopython 库。")
    print("请先通过 'pip install biopython' 命令进行安装。")
    sys.exit(1)

def write_vcf_header(output_handle, reference_id, sample_names):
    """
    向输出文件写入 VCF 头部信息。
    """
    # 获取当前日期
    file_date = datetime.now().strftime("%Y%m%d")
    
    # 写入 VCF 版本和其他元信息
    output_handle.write("##fileformat=VCFv4.1\n")
    output_handle.write("##fileDate={}\n".format(file_date))
    output_handle.write("##source=fasta_to_vcf.py\n")
    
    # 写入参考序列的 contig 信息
    output_handle.write("##contig=<ID={}>\n".format(reference_id))
    
    # 写入 INFO 和 FILTER 字段的定义 (简化版)
    output_handle.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">\n')
    output_handle.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n')
    output_handle.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    
    # 写入 FORMAT 字段的定义，这里我们只生成基因型(GT)
    output_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    
    # 写入列标题
    header_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + sample_names
    output_handle.write("\t".join(header_columns) + "\n")

def fasta_to_vcf(fasta_path, vcf_path, ref_id):
    """
    执行从 FASTA 到 VCF 转换的核心函数。
    """
    print("开始处理 FASTA 文件: {}".format(fasta_path))
    
    # 1. 读取所有 FASTA 记录
    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            print("错误: FASTA 文件为空或格式不正确。")
            return
    except FileNotFoundError:
        print("错误: 找不到输入文件 '{}'".format(fasta_path))
        return
    
    # 2. 查找参考序列和其他样本序列
    reference_record = None
    sample_records = []
    
    for record in records:
        if record.id == ref_id:
            reference_record = record
        else:
            sample_records.append(record)
            
    if reference_record is None:
        print("错误: 在 FASTA 文件中找不到指定的参考序列 ID '{}'".format(ref_id))
        print("可用的序列 ID 有: {}".format(", ".join([rec.id for rec in records])))
        return
        
    print("找到参考序列: {}".format(reference_record.id))
    print("找到 {} 个样本序列进行比较。".format(len(sample_records)))

    # 3. 检查所有序列长度是否一致
    alignment_length = len(reference_record.seq)
    for record in sample_records:
        if len(record.seq) != alignment_length:
            print("警告: 序列 '{}' 的长度 ({}) 与参考序列 ({}) 不一致。请确保输入是比对好的文件。".format(
                record.id, len(record.seq), alignment_length
            ))
            # 可以选择在这里退出或继续处理
            # sys.exit(1)

    # 4. 打开输出文件并写入 VCF 头部
    with open(vcf_path, 'w') as f_out:
        # VCF 文件中的样本列应包含所有序列，包括参考序列
        all_sample_ids = [rec.id for rec in [reference_record] + sample_records]
        write_vcf_header(f_out, reference_record.id, all_sample_ids)
        
        ref_seq = str(reference_record.seq).upper()
        
        # 5. 遍历比对的每一个位置
        for i in range(alignment_length):
            pos = i + 1
            ref_base = ref_seq[i]
            
            # 如果参考位点是 gap 或 'N'，则跳过
            if ref_base in ('-', 'N'):
                continue
                
            # 获取该位置上所有样本的碱基，并找出所有不同于参考的碱基(ALT)
            position_bases = {str(rec.seq[i]).upper() for rec in sample_records}
            alt_alleles = sorted(list(position_bases - {ref_base, '-', 'N'}))
            
            # 如果存在变异，则写入 VCF 记录
            if alt_alleles:
                alt_string = ",".join(alt_alleles)
                
                # 准备 INFO 字段
                allele_counts = []
                total_alleles = 0

                # 生成每个样本的基因型
                genotypes = []
                for sample_rec in [reference_record] + sample_records:
                    sample_base = str(sample_rec.seq[i]).upper()
                    
                    if sample_base == ref_base:
                        # 与参考基因型相同
                        genotype = '0/0'
                        total_alleles += 2
                    elif sample_base in alt_alleles:
                        # 是一个 ALT 等位基因
                        alt_index = alt_alleles.index(sample_base) + 1
                        genotype = '{0}/{0}'.format(alt_index)
                        # 确保 allele_counts 列表足够长
                        while len(allele_counts) < alt_index:
                            allele_counts.append(0)
                        allele_counts[alt_index-1] += 2
                        total_alleles += 2
                    else:
                        # 缺失数据 (gap or N)
                        genotype = './.'
                    genotypes.append(genotype)

                # 构建 INFO 字段字符串
                ac_str = "AC=" + ",".join(map(str, allele_counts)) if allele_counts else "."
                an_str = "AN=" + str(total_alleles)
                info_str = "{};{}".format(ac_str, an_str)

                # 准备 VCF 行的其他字段
                chrom = reference_record.id
                vcf_id = "."
                qual = "." # FASTA 文件没有质量信息
                filter_status = "PASS"
                vcf_format = "GT"
                
                # 组合成一行并写入文件
                vcf_line = [
                    chrom, str(pos), vcf_id, ref_base, alt_string,
                    qual, filter_status, info_str, vcf_format
                ] + genotypes
                f_out.write("\t".join(vcf_line) + "\n")

    print("\n转换完成！ VCF 文件已保存至: {}".format(vcf_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="将 FASTA 格式的多序列比对文件转换为 VCF 格式。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="输入的 FASTA 文件路径。"
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="输出的 VCF 文件路径。"
    )
    parser.add_argument(
        "-r", "--ref", 
        required=True, 
        help="在 FASTA 文件中用作参考序列的序列 ID (例如 >sequence1 中的 \"sequence1\")。"
    )
    
    args = parser.parse_args()
    
    fasta_to_vcf(args.input, args.output, args.ref)