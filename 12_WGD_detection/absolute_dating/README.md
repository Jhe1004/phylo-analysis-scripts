# 12_WGD_detection/absolute_dating

## 功能概述
这是一条通过基因树定年和 duplication age 分布来识别潜在 WGD 峰的流程。

## 子目录说明
- `0_get_gene_family`：整理基因家族
- `1_mafft`：批量比对
- `2_delmissing`：过滤缺失位点
- `3_raxml`：构建基因树
- `4_reroot_gene_family_trees`：重新定根
- `5_extract_time`：提取 duplication age
- `6_generate_treePL_file`：生成 treePL 配置
- `7_visualize_results`：绘图和结果汇总

## 使用建议
请优先进入具体子目录阅读对应的 `README.md`，再修改脚本顶部配置并运行。
