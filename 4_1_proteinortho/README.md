# 4_1_proteinortho

## 功能概述
用于基于蛋白序列运行 ProteinOrtho，识别多物种之间的同源基因簇。脚本会先准备输入蛋白，再调用 ProteinOrtho 和比对程序，输出可供下游抽取序列和系统发育分析使用的 ortholog 表格。

## 正式入口
- `run_proteinortho.py`

## 输入与输出
### 输入
- `input/` 中各物种蛋白序列文件
### 输出
- `output/` 中的 ProteinOrtho 结果表

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
- 优先使用 `diamond`，找不到时回退到 `blastp`。
