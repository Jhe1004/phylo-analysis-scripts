# 20_Chloroplast/9_复杂基因结构建树

## 功能概述
这一部分面向复杂叶绿体基因结构的系统发育分析，包含从片段拆分、同源识别、比对过滤到 IQ-TREE / ExaBayes 建树的一整套步骤。

## 子目录说明
- `0_把gb分割成小片段`：从 gb 文件中拆分片段
- `1_proteinortho`：识别同源关系
- `2_split_tsv_by_type`：拆分 orthology 结果
- `3_parse_ortho_dna_groups`：整理 DNA 同源组
- `4_mafft`：比对并可选拼接矩阵
- `5_delmissing`：过滤缺失位点
- `6_iqtree`：运行 IQ-TREE
- `7_exabayes`：运行 ExaBayes

## 使用建议
请优先进入具体子目录阅读对应的 `README.md`，再修改脚本顶部配置并运行。
