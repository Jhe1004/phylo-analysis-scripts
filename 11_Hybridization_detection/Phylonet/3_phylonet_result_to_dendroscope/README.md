# Phylonet Step 3

## 功能概述
从 PhyloNet 结果文本中抽取 Dendroscope 可视化树和总 log probability，方便后续挑选最佳网络并展示。

## 正式入口
- `get_tree_and_probability_from_phylonet.py`

## 输入与输出
### 输入
- `input/*.txt`
### 输出
- `output/tree/`
- `output/probability/result.txt`

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
