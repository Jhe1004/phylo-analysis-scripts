# 1_trinity

## 功能概述
用于批量完成 Trinity 转录组组装，并在组装结束后自动提取每个样本的最长转录本。这个目录已经整理成一个主调度脚本，内部顺序调用 helper 脚本。

## 正式入口
- `trinity_pipeline.py`

## 输入与输出
### 输入
- `input/` 中的 clean reads
### 输出
- `output/assemblies/` 中的 Trinity 组装结果
- `output/longest_isoforms/` 中的最长转录本 FASTA

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
- 依赖 helper 脚本位于 `dependencies/scripts/`。
- 适合作为 Trimmomatic 之后的第二步。
