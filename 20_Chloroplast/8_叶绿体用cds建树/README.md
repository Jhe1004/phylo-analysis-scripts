# 20_Chloroplast/8_叶绿体用cds建树

## 功能概述
用于从叶绿体 GenBank 文件中抽取同源 CDS，整理成后续比对和建树可直接使用的矩阵。

## 正式入口
- `extract_orthologs.py`

## 输入与输出
### 输入
- `input/` 中的参考与目标 GenBank 文件
### 输出
- `output/` 中的同源 CDS FASTA

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
