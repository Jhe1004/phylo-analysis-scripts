# 3_TransDecoder-master

## 功能概述
用于对转录本序列进行 ORF 预测，提取潜在编码区并生成 `.pep` 和 `.cds` 等结果文件。脚本会按样本建立独立工作目录，方便管理中间文件和最终结果。

## 正式入口
- `transdecoder.py`

## 输入与输出
### 输入
- `input/` 中待预测的转录本 FASTA
### 输出
- `output/work/样本名/` 中的 TransDecoder 运行结果

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
- 会自动查找 `TransDecoder.LongOrfs` 和 `TransDecoder.Predict`。
- 上游原始工具文件保留在 `dependencies/upstream/`。
