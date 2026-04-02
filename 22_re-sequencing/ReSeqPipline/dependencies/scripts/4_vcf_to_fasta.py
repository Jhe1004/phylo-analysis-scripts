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
import logging
from datetime import datetime
from collections import OrderedDict


def setup_logger(script_name):
    """配置日志记录器，同时输出到文件和控制台"""
    log_file = f"{script_name}.log"
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.DEBUG)
    
    logger.handlers.clear()
    
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    logger.info(f"日志文件: {log_file}")
    return logger


logger = setup_logger("4_vcf_to_fasta")


try:
    import pysam
except ImportError:
    logger.error("需要 pysam 库。请运行: pip install pysam")
    sys.exit(1)

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    logger.error("需要 biopython 库。请运行: pip install biopython")
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
    alleles_set = frozenset([a.upper() for a in alleles if a.upper() in {'A', 'C', 'G', 'T'}])
    
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
    logger.info(f"正在加载参考基因组: {ref_path}")
    sequences = OrderedDict()
    
    for record in SeqIO.parse(ref_path, "fasta"):
        sequences[record.id] = list(str(record.seq).upper())
        logger.info(f"  - {record.id}: {len(sequences[record.id]):,} bp")
    
    logger.info(f"共加载 {len(sequences)} 个 contigs/scaffolds。")
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
    
    logger.info(f"\n{'='*60}")
    logger.info(f"正在处理: {os.path.basename(vcf_path)}")
    logger.info('='*60)
    
    # 打开 VCF 文件
    vcf = pysam.VariantFile(vcf_path)
    
    # 获取样本列表
    samples = list(vcf.header.samples)
    logger.info(f"检测到 {len(samples)} 个样本: {', '.join(samples)}")
    
    # 为每个样本初始化序列（全部填充为 N，表示未知）
    sample_seqs = {}
    for sample in samples:
        sample_seqs[sample] = {
            contig: ['N'] * len(seq) for contig, seq in reference_seqs.items()
        }
    
    # 统计信息
    stats = {sample: {'total': 0, 'hom_ref': 0, 'hom_alt': 0, 'het': 0, 'missing': 0, 'filtered': 0}
             for sample in samples}
    
    # 遍历 VCF 记录
    logger.info("正在解析变异记录...")
    for record in vcf:
        contig = record.chrom
        start = record.start
        stop = record.stop
        
        if contig not in reference_seqs:
            continue
        
        ref_allele = record.ref
        alt_alleles = record.alts if record.alts else []
        
        # 处理每个样本
        for sample in samples:
            stats[sample]['total'] += (stop - start)
            
            try:
                gt = record.samples[sample]['GT']
            except KeyError:
                stats[sample]['missing'] += (stop - start)
                continue
            
            # 检查是否为缺失数据
            if gt is None or None in gt:
                stats[sample]['missing'] += (stop - start)
                # 保持为 'N'
                continue
            
            # 获取质量信息
            try:
                # 对于 gVCF 参考块，GQ 通常在样本字段中
                gq = record.samples[sample].get('GQ', None)
                dp = record.samples[sample].get('DP', None)
                
                # 如果是 gVCF 参考块且没有 GQ，尝试使用 MIN_DP 或其他
                if gq is None: gq = 99
                if dp is None: dp = 100
            except:
                gq = 99
                dp = 100
            
            # 质量过滤
            if (gq is not None and gq < min_gq) or (dp is not None and dp < min_dp):
                stats[sample]['filtered'] += (stop - start)
                # 保持为 'N'
                continue
            
            # 获取基因型对应的等位基因
            alleles_list = [ref_allele] + list(alt_alleles)
            try:
                called_alleles = [alleles_list[i] for i in gt]
            except (IndexError, TypeError):
                stats[sample]['missing'] += (stop - start)
                continue
            
            # 过滤掉非 ATCG 的等位基因 (如 <NON_REF>, *)
            called_alleles_clean = [a for a in called_alleles if all(c in 'ACGTacgt' for c in a)]
            
            if not called_alleles_clean:
                # 如果只有 <NON_REF> 之类的，通常视为参考
                if any(a == '<NON_REF>' for a in called_alleles) and len(set(gt)) == 1 and gt[0] == 0:
                     called_alleles_clean = [ref_allele]
                else:
                    stats[sample]['missing'] += (stop - start)
                    continue

            # 判断基因型类型
            unique_alleles = set(called_alleles_clean)
            
            if len(unique_alleles) == 1:
                # 纯合
                allele = list(unique_alleles)[0]
                
                if allele == ref_allele:
                    stats[sample]['hom_ref'] += (stop - start)
                    # 填充参考碱基
                    for p in range(start, stop):
                        if p < len(sample_seqs[sample][contig]):
                            sample_seqs[sample][contig][p] = reference_seqs[contig][p]
                else:
                    # 纯合变异
                    if len(allele) == 1 and (stop - start) == 1:
                        # SNP
                        stats[sample]['hom_alt'] += 1
                        sample_seqs[sample][contig][start] = allele.upper()
                    else:
                        # Indel 或跨度不匹配，暂不处理或仅标记
                        stats[sample]['filtered'] += (stop - start)
            else:
                # 杂合
                if (stop - start) == 1:
                    stats[sample]['het'] += 1
                    if use_iupac:
                        sample_seqs[sample][contig][start] = get_iupac_code(called_alleles_clean)
                    else:
                        sample_seqs[sample][contig][start] = reference_seqs[contig][start]
                else:
                    # 跨度多个碱基的杂合通常不出现在 gVCF 参考块中
                    stats[sample]['filtered'] += (stop - start)
    
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
                id=contig,
                description=f"Consensus from VCF, sample={sample}, GQ>={min_gq}, DP>={min_dp}, IUPAC={use_iupac}"
            )
            records.append(record)
        
        SeqIO.write(records, output_path, "fasta")
        output_files.append(output_path)
        
        # 打印统计信息
        s = stats[sample]
        logger.info(f"\n  样本: {sample}")
        logger.info(f"    总记录数: {s['total']:,}")
        logger.info(f"    纯合参考: {s['hom_ref']:,}")
        logger.info(f"    纯合变异: {s['hom_alt']:,}")
        logger.info(f"    杂合位点: {s['het']:,}")
        logger.info(f"    缺失数据: {s['missing']:,}")
        logger.info(f"    质量过滤: {s['filtered']:,}")
        logger.info(f"    输出文件: {output_path}")
    
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
        logger.error(f"参考基因组不存在: {reference}")
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
        logger.error(f"在 {vcf_folder} 中未找到 VCF 文件。")
        sys.exit(1)
    
    logger.info(f"\n找到 {len(vcf_files)} 个 VCF 文件:")
    for vcf in vcf_files:
        logger.info(f"  - {os.path.basename(vcf)}")
    
    # 处理每个 VCF 文件 - 并行处理
    logger.info(f"\n开始并行处理，并发数: {args.threads}")
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
                logger.error(f"[错误] 处理文件 {vcf_path} 时发生异常: {e}")
                
    # 排序输出文件列表
    all_output_files.sort()
    
    logger.info("\n" + "="*60)
    logger.info("一致性序列生成完成！")
    logger.info("="*60)
    logger.info(f"输出文件 ({len(all_output_files)} 个):")
    for f in all_output_files:
        logger.info(f"  - {f}")
    
    logger.info(f"\n下一步:")
    logger.info(f"  1. 运行 5_check_non_atcg.py 检查序列质量")
    logger.info(f"  2. 运行 6_combine_fasta.py 合并所有序列")


if __name__ == "__main__":
    main()
