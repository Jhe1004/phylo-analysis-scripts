# phypart Step 0

## 功能概述
从直系同源结果中整理 gene family，为基于 PhyParts 的 WGD 检测准备输入。

## 正式入口
- `get_gene_family.py`

## 输入与输出
### 输入
- `input/` 中的同源结果和序列文件
### 输出
- `output/` 中的基因家族序列

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
