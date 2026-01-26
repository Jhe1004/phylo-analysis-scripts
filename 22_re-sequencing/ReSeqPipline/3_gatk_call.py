#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用 GATK HaplotypeCaller 进行变异检测。

本脚本执行以下步骤：
1. 对每个样本的 BAM 文件运行 GATK HaplotypeCaller，生成 gVCF 文件
2. (可选) 使用 CombineGVCFs + GenotypeGVCFs 合并多个样本
3. 或者直接生成单样本 VCF

对于系统发育分析，我们使用 --emit-ref-confidence GVCF 模式，
这样可以区分"没有变异"和"没有数据覆盖"的位点。

软件依赖:
- gatk (版本 4.x)
"""

import os
import argparse
import subprocess
import sys
import glob


def run_command(command, description, exit_on_fail=True):
    """运行命令并打印状态"""
    print(f"\n  [子步骤] {description}")
    print(f"  [命令] {command}")
    
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print(f"  [错误] 命令执行失败，返回码: {result.returncode}")
        if exit_on_fail:
            sys.exit(1)
        return False
    return True


def call_variants_single_sample(bam_path, reference, output_folder, emit_gvcf=True):
    """对单个样本进行变异检测"""
    
    sample_name = os.path.basename(bam_path).replace('.bam', '')
    
    print(f"\n{'='*60}")
    print(f"[样本] {sample_name}")
    print('='*60)
    
    if emit_gvcf:
        output_file = os.path.join(output_folder, f"{sample_name}.g.vcf.gz")
        erc_option = "--emit-ref-confidence GVCF"
        desc = "GATK HaplotypeCaller (生成 gVCF)"
    else:
        output_file = os.path.join(output_folder, f"{sample_name}.vcf.gz")
        erc_option = ""
        desc = "GATK HaplotypeCaller (生成 VCF)"
    
    # GATK HaplotypeCaller 命令
    # --native-pair-hmm-threads: 使用多线程加速（需要 GATK 4.1+）
    hc_cmd = (
        f"gatk HaplotypeCaller "
        f"-R {reference} "
        f"-I {bam_path} "
        f"-O {output_file} "
        f"{erc_option} "
        f"--native-pair-hmm-threads 4"
    )
    
    run_command(hc_cmd, desc)
    
    print(f"[完成] 输出: {output_file}")
    return output_file


def combine_gvcfs(gvcf_files, reference, output_folder):
    """合并多个 gVCF 文件"""
    
    print(f"\n{'='*60}")
    print("[步骤] 合并所有样本的 gVCF 文件")
    print('='*60)
    
    combined_gvcf = os.path.join(output_folder, "cohort.g.vcf.gz")
    
    # 构建 -V 参数列表
    v_args = " ".join([f"-V {gvcf}" for gvcf in gvcf_files])
    
    combine_cmd = (
        f"gatk CombineGVCFs "
        f"-R {reference} "
        f"{v_args} "
        f"-O {combined_gvcf}"
    )
    
    run_command(combine_cmd, "GATK CombineGVCFs")
    
    return combined_gvcf


def genotype_gvcfs(combined_gvcf, reference, output_folder):
    """对合并的 gVCF 进行基因型分析"""
    
    print(f"\n{'='*60}")
    print("[步骤] 进行联合基因型分析 (Joint Genotyping)")
    print('='*60)
    
    final_vcf = os.path.join(output_folder, "all_samples.vcf.gz")
    
    genotype_cmd = (
        f"gatk GenotypeGVCFs "
        f"-R {reference} "
        f"-V {combined_gvcf} "
        f"-O {final_vcf} "
        f"--include-non-variant-sites"  # 重要：包含非变异位点，用于系统发育
    )
    
    run_command(genotype_cmd, "GATK GenotypeGVCFs")
    
    print(f"\n[完成] 联合 VCF 文件: {final_vcf}")
    return final_vcf


def main():
    parser = argparse.ArgumentParser(
        description="使用 GATK HaplotypeCaller 对 BAM 文件进行变异检测。",
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
        '-b', '--bam-folder',
        default=".",
        metavar='DIR',
        help="包含 BAM 文件的文件夹 (默认: 当前目录)"
    )
    optional.add_argument(
        '-o', '--output',
        default="./vcf_output",
        metavar='DIR',
        help="VCF 输出文件夹 (默认: ./vcf_output)"
    )
    optional.add_argument(
        '--joint-calling',
        action='store_true',
        help="进行联合变异检测 (先生成 gVCF，再合并)"
    )
    optional.add_argument(
        '--single-vcf',
        action='store_true',
        help="为每个样本生成独立的 VCF（而非 gVCF）"
    )
    
    args = parser.parse_args()
    
    # 处理路径
    bam_folder = os.path.abspath(args.bam_folder)
    output_folder = os.path.abspath(args.output)
    reference = os.path.abspath(args.reference)
    
    # 检查参考基因组及其索引
    if not os.path.exists(reference):
        print(f"错误: 参考基因组不存在: {reference}")
        sys.exit(1)
    
    dict_file = reference.rsplit('.', 1)[0] + ".dict"
    if not os.path.exists(dict_file):
        print(f"错误: 缺少 GATK 序列字典文件: {dict_file}")
        print("请先运行 1_index_reference.py 创建索引。")
        sys.exit(1)
    
    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)
    
    # 查找 BAM 文件
    bam_pattern = os.path.join(bam_folder, "*.bam")
    bam_files = sorted(glob.glob(bam_pattern))
    
    # 排除中间文件
    bam_files = [b for b in bam_files if not any(x in b for x in ['.raw.', '.sorted.'])]
    
    if not bam_files:
        print(f"错误: 在 {bam_folder} 中未找到 BAM 文件。")
        sys.exit(1)
    
    print(f"\n找到 {len(bam_files)} 个 BAM 文件:")
    for bam in bam_files:
        print(f"  - {os.path.basename(bam)}")
    
    # 决定是生成 gVCF 还是普通 VCF
    emit_gvcf = args.joint_calling or (not args.single_vcf)
    
    # 对每个样本进行变异检测
    output_files = []
    for bam_path in bam_files:
        vcf_path = call_variants_single_sample(
            bam_path, 
            reference, 
            output_folder,
            emit_gvcf=emit_gvcf
        )
        output_files.append(vcf_path)
    
    # 如果是联合检测模式，合并 gVCF
    if args.joint_calling and len(output_files) > 1:
        combined_gvcf = combine_gvcfs(output_files, reference, output_folder)
        final_vcf = genotype_gvcfs(combined_gvcf, reference, output_folder)
        
        print("\n" + "="*60)
        print("变异检测完成！")
        print("="*60)
        print(f"联合 VCF 文件: {final_vcf}")
    else:
        print("\n" + "="*60)
        print("变异检测完成！")
        print("="*60)
        print(f"输出文件 ({len(output_files)} 个):")
        for f in output_files:
            print(f"  - {f}")
    
    print("\n下一步: 运行 4_vcf_to_fasta.py 从 VCF 生成一致性序列。")


if __name__ == "__main__":
    main()
