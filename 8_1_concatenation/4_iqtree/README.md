# 8_1_concatenation/4_iqtree

## 功能概述
用于对拼接矩阵或多个矩阵运行 IQ-TREE。脚本支持自动扫描输入、可选自动拼接小矩阵，以及联合矩阵的 partition 生成。

## 正式入口
- `run_iqtree_batch_simple.py`

## 输入与输出
### 输入
- `input/` 中的 FASTA / NEXUS 矩阵
### 输出
- `output/` 中的 IQ-TREE 结果

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
- 支持 `INPUT_EXTENSIONS`、`PARTITION_FILE` 和 `AUTO_CONCAT_SMALL_MATRICES` 等配置。
