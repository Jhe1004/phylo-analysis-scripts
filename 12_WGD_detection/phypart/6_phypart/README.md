# phypart Step 6

## 功能概述
调用 PhyParts 对 gene tree 与 species tree 的冲突关系进行统计，从而辅助识别 WGD 和基因重复信号。

## 正式入口
- `run_phypart.py`

## 输入与输出
### 输入
- `input/` 中的 gene trees 和 species tree
### 输出
- `output/` 中的 PhyParts 统计结果

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
- 当前需要你自行准备可用的 PhyParts jar。
