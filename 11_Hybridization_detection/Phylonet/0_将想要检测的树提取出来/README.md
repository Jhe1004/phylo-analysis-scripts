# Phylonet Step 0

## 功能概述
从总基因树集合中提取包含目标物种组合的子树，并统一用指定外类群重新定根。这个步骤的输出是后续生成 PhyloNet 输入文件的基础。

## 正式入口
- `extract_subtrees.py`

## 输入与输出
### 输入
- `input/result.tree`
- `input/list.txt`
### 输出
- `output/out.trees`

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
