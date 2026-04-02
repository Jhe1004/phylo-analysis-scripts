# 8_genetree_speciestree_conflict/0_extract_subtree

## 功能概述
用于从大批基因树中提取关注类群的子树，作为后续冲突分析和 densitree 可视化的输入。

## 正式入口
- `extract_subtrees.py`

## 输入与输出
### 输入
- `input/` 中的基因树文件和目标物种列表
### 输出
- `output/` 中提取后的子树集合

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
