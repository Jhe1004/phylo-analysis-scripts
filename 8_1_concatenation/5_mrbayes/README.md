# 8_1_concatenation/5_mrbayes

## 功能概述
用于对拼接矩阵运行 MrBayes 贝叶斯建树。脚本支持直接读取联合矩阵、调用分区文件，并管理 MrBayes 运行参数。

## 正式入口
- `run_mrbayes_workflow.py`

## 输入与输出
### 输入
- `input/` 中的拼接矩阵与分区文件
### 输出
- `output/` 中的 MrBayes 结果和日志

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
