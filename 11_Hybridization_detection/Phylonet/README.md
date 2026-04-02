# 11_Hybridization_detection/Phylonet

## 功能概述
这一部分用于使用 PhyloNet 进行杂交网络推断。流程被拆成 4 个独立步骤，方便逐步检查中间结果。

## 子目录说明
- `0_将想要检测的树提取出来`：从总基因树中提取目标子树
- `1_Get_phylonet_input`：生成 `.nex` 输入文件
- `2_Run_phylonet`：并行运行 PhyloNet
- `3_phylonet_result_to_dendroscope`：整理网络结果并提取可视化树

## 使用建议
请优先进入具体子目录阅读对应的 `README.md`，再修改脚本顶部配置并运行。
