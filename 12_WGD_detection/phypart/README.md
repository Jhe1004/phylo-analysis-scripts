# 12_WGD_detection/phypart

## 功能概述
这是一条通过 gene tree 与 species tree 冲突统计来辅助定位 WGD 的流程。

## 子目录说明
- `0_get_gene_family`：整理基因家族
- `1_mafft`：批量比对
- `2_del_indel`：过滤缺失位点
- `3_raxml`：构建 gene tree
- `4_reroot_gene_family_trees`：重新定根
- `5_rename`：整理叶名格式
- `6_phypart`：运行 PhyParts

## 使用建议
请优先进入具体子目录阅读对应的 `README.md`，再修改脚本顶部配置并运行。
