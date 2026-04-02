# 4_2_get_seq_from_proteinortho

## 功能概述
用于根据 ProteinOrtho 的结果表，从对应的 `.pep` 和 `.cds` 文件中抽取同源序列，整理成后续比对和建树所需的基因矩阵。

## 正式入口
- `get_seq_from_proteinortho.py`

## 输入与输出
### 输入
- `input/proteinortho.tsv`
- `input/sequences/` 中的 `.pep` 与 `.cds` 文件
### 输出
- `output/` 中按基因拆分好的蛋白和核酸序列文件

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
- 保留了直系同源筛选和结果汇总逻辑。
