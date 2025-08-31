# -*- coding: utf-8 -*-
import os
import sys
import csv

# 尝试导入 dendropy，如果失败则提示安装
try:
    import dendropy
except ImportError:
    print("错误：需要 'dendropy' 库。请使用 'pip install dendropy' 命令安装它。")
    sys.exit(1)

# --- 用户配置区域 ---

# 1. 包含 Newick 树文件的文件夹路径
TREE_DIR = './'

# 2. 树文件名的后缀 (例如，RAxML 输出的 RAxML_bipartitions.result)
TREE_SUFFIX = '.treefile'

# 3. 关键物种名称列表
SPECIES_A_LIST = [
    "Clematis_repens.fasta.transdecoder.pep",
    "Clematis_otophora.fasta.transdecoder.pep"  # 示例，可替换为实际名称
]

SPECIES_B_LIST = [
    "Clematis_songorica_kashi.fasta.transdecoder.pep","Clematis_songorica_aleotai.fasta.transdecoder.pep","Clematis_songorica_hami.fasta.transdecoder.pep","Clematis_songorica_wulumuqi.fasta.transdecoder.pep"    # 示例，可替换为实际名称
]

# 4. 需要统计亲缘关系的物种列表 (SOI)
SPECIES_OF_INTEREST = [
    "Clematis_glauca_liancheng.fasta.transdecoder.pep","Clematis_intricata_2019052001.fasta.transdecoder.pep","Clematis_intricata_weijing.fasta.transdecoder.pep","Clematis_akebioides_wenchuan.fasta.transdecoder.pep","Clematis_akebioides_ganzi.fasta.transdecoder.pep","Clematis_tangutica_yongchang.fasta.transdecoder.pep","Clematis_tangutica_maerkang.fasta.transdecoder.pep","Clematis_tenuifolia.fasta.transdecoder.pep"
    # 在这里添加更多物种名称
]

# 5. 节点支持率阈值 (Bootstrap Value)
SUPPORT_THRESHOLD = 50

# 6. 输出 CSV 文件的名称
OUTPUT_CSV_FILE = 'relationship_results_dendropy.csv'

# --- 脚本主逻辑 ---

def analyze_relationships(tree_dir, tree_suffix, species_a_list, species_b_list, soi_list, support_threshold, output_csv_file):
    """
    使用 dendropy 分析亲缘关系并将结果保存到 CSV 文件。
    """
    results = []
    total_files_processed = 0
    total_trees_analyzed = 0

    print(f"开始分析目录 '{tree_dir}' 中以 '{tree_suffix}' 结尾的树文件...")
    print(f"支持率阈值设置为: > {support_threshold}")
    print(f"物种 A 列表: {', '.join(species_a_list)}")
    print(f"物种 B 列表: {', '.join(species_b_list)}")
    print(f"待分析物种 (SOI): {', '.join(soi_list)}")
    print("-" * 30)

    if not os.path.isdir(tree_dir):
        print(f"错误: 目录 '{tree_dir}' 不存在或不是一个有效的目录。")
        sys.exit(1)

    for filename in os.listdir(tree_dir):
        if filename.endswith(tree_suffix):
            filepath = os.path.join(tree_dir, filename)
            total_files_processed += 1
            print(f"\n正在处理文件: {filename}")

            try:
                # 使用 dendropy 加载树文件
                tree = dendropy.Tree.get(path=filepath, schema="newick")
            except Exception as e:
                print(f"错误: 使用 dendropy 解析树文件 '{filepath}' 失败。错误信息: {e}")
                continue

            total_trees_analyzed += 1
            
            # 获取树中存在的所有物种标签，用于快速检查
            taxon_labels_in_tree = {t.label for t in tree.taxon_namespace}

            for a_name in species_a_list:
                for b_name in species_b_list:
                    for soi_name in soi_list:
                        
                        # 首先检查所有三个物种是否存在于当前树中
                        required_labels = {a_name, b_name, soi_name}
                        if not required_labels.issubset(taxon_labels_in_tree):
                            continue

                        relationship = 'error' # 默认值
                        try:
                            # 核心逻辑：比较三个最近共同祖先节点
                            mrca_sab = tree.mrca(taxon_labels=[soi_name, a_name, b_name])
                            mrca_sa = tree.mrca(taxon_labels=[soi_name, a_name])
                            mrca_sb = tree.mrca(taxon_labels=[soi_name, b_name])

                            if mrca_sa is None or mrca_sb is None or mrca_sab is None:
                                relationship = 'error' # 无法找到祖先
                            elif mrca_sa == mrca_sb:
                                # (A, B, SOI) 形成多分叉，关系不确定
                                relationship = 'undetermined'
                            elif mrca_sab == mrca_sb:
                                # (SOI, B)聚为一枝，A是其外类群，说明SOI离B更近
                                # 从 node.label 读取支持率
                                support_str = mrca_sb.label
                                if support_str is not None and float(support_str) > support_threshold:
                                    relationship = 'closer_to_B'
                                else:
                                    relationship = 'undetermined' # 支持率不足
                            elif mrca_sab == mrca_sa:
                                # (SOI, A)聚为一枝，B是其外类群，说明SOI离A更近
                                # 从 node.label 读取支持率
                                support_str = mrca_sa.label
                                if support_str is not None and float(support_str) > support_threshold:
                                    relationship = 'closer_to_A'
                                else:
                                    relationship = 'undetermined' # 支持率不足
                            else:
                                # 罕见的拓扑结构，无法判断
                                relationship = 'undetermined'
                        
                        except Exception as e_mrca:
                            print(f"  错误: 为 SOI '{soi_name}' 与 A '{a_name}' 和 B '{b_name}' 分析时出错: {e_mrca}")
                            relationship = 'error'

                        # 记录结果
                        result_entry = {
                            'SOI': soi_name,
                            'A': a_name,
                            'B': b_name,
                            'tree': filename,
                            'relationship': relationship
                        }
                        results.append(result_entry)

    # --- 汇总统计 (与之前版本相同) ---
    summary = {}
    for soi in soi_list:
        for a in species_a_list:
            for b in species_b_list:
                key = (soi, a, b)
                summary[key] = {
                    'trees_processed': 0,
                    'closer_to_A': 0,
                    'closer_to_B': 0,
                    'undetermined': 0,
                    'errors': 0
                }

    for res in results:
        key = (res['SOI'], res['A'], res['B'])
        # 确保在汇总前所有物种都存在
        if key not in summary:
            continue
        summary[key]['trees_processed'] += 1
        if res['relationship'] == 'closer_to_A':
            summary[key]['closer_to_A'] += 1
        elif res['relationship'] == 'closer_to_B':
            summary[key]['closer_to_B'] += 1
        elif res['relationship'] == 'undetermined':
            summary[key]['undetermined'] += 1
        elif res['relationship'] == 'error':
            summary[key]['errors'] += 1

    # --- 打印和输出结果 (与之前版本相同) ---
    print("\n" + "=" * 30)
    print("分析完成。结果摘要:")
    print(f"总共扫描到符合后缀的文件数: {total_files_processed}")
    print(f"成功加载并分析的树文件数: {total_trees_analyzed}")
    print("-" * 30)

    for key, counts in summary.items():
        soi, a, b = key
        print(f"SOI: {soi}, A: {a}, B: {b}")
        print(f"  在 A, B, SOI 均存在的树中的分析次数: {counts['trees_processed']}")
        print(f"  -> 更接近 A '{a}' (支持率 > {support_threshold}) 的次数: {counts['closer_to_A']}")
        print(f"  -> 更接近 B '{b}' (支持率 > {support_threshold}) 的次数: {counts['closer_to_B']}")
        print(f"  -> 关系不确定/等距/支持率不足的次数: {counts['undetermined']}")
        print(f"  处理中遇到错误的次数: {counts['errors']}")
        print("-" * 20)

    print(f"\n正在将结果写入 CSV 文件: {output_csv_file} ...")
    try:
        with open(output_csv_file, 'w', newline='', encoding='utf-8') as csvfile:
            header = [
                'SOI',
                'A',
                'B',
                'Trees Analyzed (A, B, SOI present)',
                f'Closer to A (Support > {support_threshold})',
                f'Closer to B (Support > {support_threshold})',
                'Undetermined/Equidistant/Low Support',
                'Errors'
            ]
            writer = csv.writer(csvfile)
            writer.writerow(header)

            for key, counts in summary.items():
                soi, a, b = key
                row = [
                    soi, a, b,
                    counts['trees_processed'],
                    counts['closer_to_A'],
                    counts['closer_to_B'],
                    counts['undetermined'],
                    counts['errors']
                ]
                writer.writerow(row)
        print(f"结果已成功保存到 {output_csv_file}")
    except IOError as e:
        print(f"错误: 无法写入 CSV 文件 '{output_csv_file}': {e}")
    except Exception as e_csv:
        print(f"错误: 写入 CSV 时发生未知错误: {e_csv}")

# --- 运行脚本 ---
if __name__ == "__main__":
    analyze_relationships(
        tree_dir=TREE_DIR,
        tree_suffix=TREE_SUFFIX,
        species_a_list=SPECIES_A_LIST,
        species_b_list=SPECIES_B_LIST,
        soi_list=SPECIES_OF_INTEREST,
        support_threshold=SUPPORT_THRESHOLD,
        output_csv_file=OUTPUT_CSV_FILE
    )
