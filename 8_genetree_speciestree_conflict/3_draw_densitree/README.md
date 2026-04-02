# 8_genetree_speciestree_conflict/3_draw_densitree

## 功能概述
用于把大量基因树叠加为 densitree 风格可视化，并可同时叠加物种树节点支持和支持比例饼图。

## 正式入口
- `draw_densitree_workflow.py`

## 输入与输出
### 输入
- `input/` 中的 gene trees、species tree、样本名映射等文件
### 输出
- `output/` 中的 densitree 图和辅助统计结果

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
