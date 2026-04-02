# Phylonet Step 2

## 功能概述
并行调用 `java -jar PhyloNet` 运行网络推断，对多个 `.nex` 任务批量求解杂交网络。

## 正式入口
- `aln.py`

## 输入与输出
### 输入
- `input/*.nex`
### 输出
- `output/*.txt`
- `output/run_phylonet.log`

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
