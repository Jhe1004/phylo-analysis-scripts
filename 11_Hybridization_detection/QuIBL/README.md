# QuIBL

## 功能概述
用于基于基因树分支长度检测 introgression / hybridization。`run_quibl.py` 负责调用上游 QuIBL 核心脚本，`visual_quibl.py` 负责读取结果并绘制 heatmap。

## 正式入口
- `run_quibl.py / visual_quibl.py`

## 输入与输出
### 输入
- `input/gene_trees.trees` 用于主分析
- `input/concat.newick` 与 `input/quibl_results.csv` 用于可视化
### 输出
- `output/quibl_results.csv`
- `output/heatmap_with_tree.png`
- `output/run_quibl.log`

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找

## 备注
- 上游脚本位于 `dependencies/upstream/`。
