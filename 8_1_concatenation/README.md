# 8_1_concatenation

## 功能概述
这一部分用于基于拼接矩阵进行系统发育分析，从矩阵拼接、缺失位点过滤到 RAxML / IQ-TREE / MrBayes 建树都包含在内。

## 子目录说明
- `0_concat`：拼接多个单基因矩阵并生成分区文件
- `1_del_indel`：过滤拼接矩阵中的缺失位点
- `2_codon`：拆分或重组密码子矩阵
- `3_raxml`：运行 RAxML 建树
- `4_iqtree`：运行 IQ-TREE 建树
- `5_mrbayes`：运行 MrBayes 贝叶斯建树

## 使用建议
请优先进入具体子目录阅读对应的 `README.md`，再修改脚本顶部配置并运行。
