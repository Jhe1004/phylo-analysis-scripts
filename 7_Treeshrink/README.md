# 7_Treeshrink

## 功能概述
用于先为每个基因矩阵建树，再调用 TreeShrink 对异常长枝或潜在异常序列进行检测和裁剪，最后输出收缩后的比对矩阵。

## 正式入口
- `run_pipeline.py`

## 输入与输出
### 输入
- `input/alignments/` 中的基因比对文件
### 输出
- `output/gene_trees/` 中的基因树
- `output/shrunk_alignments/` 中的收缩后比对结果

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
- 建树和过滤步骤都支持并行。
- 依赖 RAxML、TreeShrink 和 `Rscript`。
