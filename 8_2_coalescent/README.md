# 8_2_coalescent

## 功能概述
这一部分用于基于 gene tree 的共祖模型分析，核心思路是先建 gene tree，再过滤，最后交给 ASTRAL 进行 species tree 推断。

## 子目录说明
- `0_raxml`：为每个基因独立建树
- `1_filter`：过滤基因矩阵和 gene tree
- `2_0_astral`：运行 ASTRAL 推断 species tree

## 使用建议
请优先进入具体子目录阅读对应的 `README.md`，再修改脚本顶部配置并运行。
