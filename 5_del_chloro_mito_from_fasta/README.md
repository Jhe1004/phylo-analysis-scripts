# 5_del_chloro_mito_from_fasta

## 功能概述
用于通过 BLAST 将叶绿体和线粒体相关序列从蛋白 FASTA 中识别并剔除，减少核基因分析中的污染或非目标序列。脚本已经改为并行版本，适合批量处理多个输入文件。

## 正式入口
- `del_chloro_mito_from_fasta.py`

## 输入与输出
### 输入
- `input/reference_gb/` 中的参考叶绿体或线粒体 GenBank 文件
- `input/sequences/` 中待过滤的蛋白 FASTA
### 输出
- `output/kept_sequences/` 中保留的序列
- `output/removed_sequences/` 中被剔除的序列
- `output/summary.tsv`

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
- 会自动查找 `makeblastdb` 和 `blastp`。
- 支持多进程并行。
