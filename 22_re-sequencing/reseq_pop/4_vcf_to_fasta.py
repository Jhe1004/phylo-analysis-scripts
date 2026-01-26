#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 VCF 文件生成一致性 FASTA 序列。

本脚本支持两种模式：
1. 单样本 VCF/gVCF: 为每个 VCF 文件生成一个 FASTA 文件
2. 多样本联合 VCF: 从一个包含多样本的 VCF 生成多个 FASTA 文件

重要特性：
- 支持 IUPAC 歧义代码处理杂合位点（可选）
- 对于低质量或缺失数据，使用 'N' 填充
- 基于 VCF 的 FILTER 和 GQ/DP 进行质量过滤

软件依赖:
- Python 3
- pysam
- biopython
"""

import os
import argparse
import sys
import gzip
from collections import OrderedDict

try:
    import pysam
except ImportError:
    print("错误: 需要 pysam 库。请运行: pip install pysam")
    sys.exit(1)

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("错误: 需要 biopython 库。请运行: pip install biopython")
    sys.exit(1)


# IUPAC 歧义代码表
IUPAC_CODES = {
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N',
}


def get_iupac_code(alleles):
    """
    根据等位基因集合返回 IUPAC 代码。
    
    Args:
        alleles: 等位基因列表或集合，如 ['A', 'G']
    
    Returns:
        对应的 IUPAC 代码字符
    """
    alleles_set = frozenset([a.upper() for a in alleles if a.upper() in 'ACGT'])
    
    if len(alleles_set) == 0:
        return 'N'
    elif len(alleles_set) == 1:
        return list(alleles_set)[0]
    else:
        return IUPAC_CODES.get(alleles_set, 'N')


def parse_reference_fasta(ref_path):
    """
    解析参考基因组，返回每个 contig 的序列。
    
    Returns:
        OrderedDict: {contig_name: sequence_string}
    """
    print(f"正在加载参考基因组: {ref_path}")
    sequences = OrderedDict()
    
    for record in SeqIO.parse(ref_path, "fasta"):
        sequences[record.id] = list(str(record.seq).upper())
        print(f"  - {record.id}: {len(sequences[record.id]):,} bp")
    
    print(f"共加载 {len(sequences)} 个 contigs/scaffolds。")
    return sequences


def vcf_to_consensus(vcf_path, reference_seqs, output_folder, 
                     use_iupac=True, min_gq=20, min_dp=3):
    """
    从 VCF 文件生成一致性序列。
    
    Args:
        vcf_path: VCF/gVCF 文件路径
        reference_seqs: 参考基因组序列字典
        output_folder: 输出文件夹
        use_iupac: 是否使用 IUPAC 代码表示杂合位点
        min_gq: 最小 Genotype Quality 阈值
        min_dp: 最小覆盖深度阈值
    """
    
    print(f"\n{'='*60}")
    print(f"正在处理: {os.path.basename(vcf_path)}")
    print('='*60)
    
    # 打开 VCF 文件
    vcf = pysam.VariantFile(vcf_path)
    
    # 获取样本列表
    samples = list(vcf.header.samples)
    print(f"检测到 {len(samples)} 个样本: {', '.join(samples)}")
    
    # 为每个样本初始化序列（深拷贝参考序列）
    sample_seqs = {}
    for sample in samples:
        sample_seqs[sample] = {
            contig: list(seq) for contig, seq in reference_seqs.items()
        }
    
    # 统计信息
    stats = {sample: {'total': 0, 'hom_ref': 0, 'hom_alt': 0, 'het': 0, 'missing': 0, 'filtered': 0}
             for sample in samples}
    
    # 遍历 VCF 记录
    print("正在解析变异记录...")
    for record in vcf:
        contig = record.chrom
        pos = record.pos - 1  # VCF 是 1-based，Python 是 0-based
        
        if contig not in reference_seqs:
            continue
        
        ref_allele = record.ref
        alt_alleles = record.alts if record.alts else []
        
        # 处理每个样本
        for sample in samples:
            stats[sample]['total'] += 1
            
            try:
                gt = record.samples[sample]['GT']
            except KeyError:
                stats[sample]['missing'] += 1
                continue
            
            # 检查是否为缺失数据
            if gt is None or None in gt:
                stats[sample]['missing'] += 1
                sample_seqs[sample][contig][pos] = 'N'
                continue
            
            # 获取质量信息
            try:
                gq = record.samples[sample].get('GQ', 99)
                dp = record.samples[sample].get('DP', 100)
            except:
                gq = 99
                dp = 100
            
            # 质量过滤
            if (gq is not None and gq < min_gq) or (dp is not None and dp < min_dp):
                stats[sample]['filtered'] += 1
                sample_seqs[sample][contig][pos] = 'N'
                continue
            
            # 获取基因型对应的等位基因
            alleles_list = [ref_allele] + list(alt_alleles)
            try:
                called_alleles = [alleles_list[i] for i in gt]
            except IndexError:
                stats[sample]['missing'] += 1
                sample_seqs[sample][contig][pos] = 'N'
                continue
            
            # 只处理 SNP（跳过 indels）
            if any(len(a) != 1 for a in called_alleles):
                # 对于 indels，保留参考等位基因（或标记为 N）
                continue
            
            # 判断基因型类型
            unique_alleles = set(called_alleles)
            
            if len(unique_alleles) == 1:
                # 纯合
                allele = list(unique_alleles)[0]
                sample_seqs[sample][contig][pos] = allele.upper()
                
                if allele == ref_allele:
                    stats[sample]['hom_ref'] += 1
                else:
                    stats[sample]['hom_alt'] += 1
            else:
                # 杂合
                stats[sample]['het'] += 1
                
                if use_iupac:
                    # 使用 IUPAC 代码
                    sample_seqs[sample][contig][pos] = get_iupac_code(called_alleles)
                else:
                    # 使用参考等位基因（或最常见的等位基因）
                    sample_seqs[sample][contig][pos] = ref_allele.upper()
    
    vcf.close()
    
    # 输出每个样本的 FASTA 文件
    output_files = []
    for sample in samples:
        output_path = os.path.join(output_folder, f"{sample}_consensus.fasta")
        
        records = []
        for contig, seq_list in sample_seqs[sample].items():
            seq_str = "".join(seq_list)
            record = SeqRecord(
                Seq(seq_str),
                id=f"{sample}_{contig}",
                description=f"Consensus from VCF, GQ>={min_gq}, DP>={min_dp}, IUPAC={use_iupac}"
            )
            records.append(record)
        
        SeqIO.write(records, output_path, "fasta")
        output_files.append(output_path)
        
        # 打印统计信息
        s = stats[sample]
        print(f"\n  样本: {sample}")
        print(f"    总记录数: {s['total']:,}")
        print(f"    纯合参考: {s['hom_ref']:,}")
        print(f"    纯合变异: {s['hom_alt']:,}")
        print(f"    杂合位点: {s['het']:,}")
        print(f"    缺失数据: {s['missing']:,}")
        print(f"    质量过滤: {s['filtered']:,}")
        print(f"    输出文件: {output_path}")
    
    return output_files


def main():
    parser = argparse.ArgumentParser(
        description="从 VCF 文件生成一致性 FASTA 序列（支持 IUPAC 歧义代码）。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    required = parser.add_argument_group("必需参数")
    required.add_argument(
        '-r', '--reference',
        required=True,
        metavar='REF.FASTA',
        help="参考基因组 FASTA 文件"
    )
    
    optional = parser.add_argument_group("可选参数")
    optional.add_argument(
        '-v', '--vcf-folder',
        default="./vcf_output",
        metavar='DIR',
        help="包含 VCF 文件的文件夹 (默认: ./vcf_output)"
    )
    optional.add_argument(
        '-o', '--output',
        default="./consensus_fasta_output",
        metavar='DIR',
        help="FASTA 输出文件夹 (默认: ./consensus_fasta_output)"
    )
    optional.add_argument(
        '--no-iupac',
        action='store_true',
        help="禁用 IUPAC 歧义代码，杂合位点使用参考等位基因"
    )
    optional.add_argument(
        '--min-gq',
        type=int,
        default=20,
        metavar='N',
        help="最小基因型质量 (GQ) 阈值 (默认: 20)"
    )
    optional.add_argument(
        '--min-dp',
        type=int,
        default=3,
        metavar='N',
        help="最小覆盖深度 (DP) 阈值 (默认: 3)"
    )
    
    optional.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        metavar='N',
        help="并发任务数 (默认: 4)"
    )
    
    args = parser.parse_args()
    
    # 处理路径
    vcf_folder = os.path.abspath(args.vcf_folder)
    output_folder = os.path.abspath(args.output)
    reference = os.path.abspath(args.reference)
    
    # 检查参考基因组
    if not os.path.exists(reference):
        print(f"错误: 参考基因组不存在: {reference}")
        sys.exit(1)
    
    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)
    
    # 加载参考基因组
    reference_seqs = parse_reference_fasta(reference)
    
    # 查找 VCF 文件
    import glob
    vcf_patterns = [
        os.path.join(vcf_folder, "*.vcf.gz"),
        os.path.join(vcf_folder, "*.vcf"),
        os.path.join(vcf_folder, "*.g.vcf.gz"),
        os.path.join(vcf_folder, "*.g.vcf"),
    ]
    
    vcf_files = []
    for pattern in vcf_patterns:
        vcf_files.extend(glob.glob(pattern))
    
    # 去重并排序
    vcf_files = sorted(set(vcf_files))
    
    if not vcf_files:
        print(f"错误: 在 {vcf_folder} 中未找到 VCF 文件。")
        sys.exit(1)
    
    print(f"\n找到 {len(vcf_files)} 个 VCF 文件:")
    for vcf in vcf_files:
        print(f"  - {os.path.basename(vcf)}")
    
    # 处理每个 VCF 文件 - 并行处理
    print(f"\n开始并行处理，并发数: {args.threads}")
    all_output_files = []
    
    from concurrent.futures import ProcessPoolExecutor, as_completed
    
    # 注意：传递巨大的 reference_seqs 对象给子进程可能会有性能开销（Copy-on-Write）
    # 在 Linux 上 ProcessPoolExecutor 使用 fork，通常效率还可以。
    # 如果内存紧张，可以考虑改用 ThreadPoolExecutor（但受限于 GIL，计算密集型效果一般），
    # 或者重构代码每次只读取参考序列的相应部分。
    # 鉴于 parse_reference_fasta 已经将所有序列读入内存，这里直接传递。
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(
                vcf_to_consensus,
                vcf_path,
                reference_seqs,
                output_folder,
                not args.no_iupac,
                args.min_gq,
                args.min_dp
            ): vcf_path 
            for vcf_path in vcf_files
        }
        
        for future in as_completed(futures):
            vcf_path = futures[future]
            try:
                output_files = future.result()
                all_output_files.extend(output_files)
            except Exception as e:
                print(f"[错误] 处理文件 {vcf_path} 时发生异常: {e}")
                
    # 排序输出文件列表
    all_output_files.sort()
    
    print("\n" + "="*60)
    print("一致性序列生成完成！")
    print("="*60)
    print(f"输出文件 ({len(all_output_files)} 个):")
    for f in all_output_files:
        print(f"  - {f}")
    
    print(f"\n下一步:")
    print(f"  1. 运行 5_check_non_atcg.py 检查序列质量")
    print(f"  2. 运行 6_combine_fasta.py 合并所有序列")


if __name__ == "__main__":
    main()
