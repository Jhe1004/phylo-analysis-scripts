# 22_re-sequencing/GetGeneFromGFF

## 功能概述
这一部分用于从 GFF 和参考基因组中提取目标基因，再进一步整理为后续分析矩阵。

## 子目录说明
- `0_filter_mRNA_annotation_from_gff`：提取目标注释
- `1_get_gene_from_fasta_in_gff`：从 FASTA 中提取基因
- `2_collection_gene`：汇总基因
- `3_del_missing`：过滤缺失基因

## 使用建议
请优先进入具体子目录阅读对应的 `README.md`，再修改脚本顶部配置并运行。
