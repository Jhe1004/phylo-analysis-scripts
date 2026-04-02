# 0_trimmomatic

## 功能概述
用于对原始双端转录组测序数据进行质控和修剪。脚本会自动匹配成对的 reads，调用 Trimmomatic 完成去接头、质量裁剪和低质量序列过滤，并输出后续 Trinity 可直接使用的 clean reads。

## 正式入口
- `trimmomatic_pe.py`

## 输入与输出
### 输入
- `input/` 中的原始双端 FASTQ 文件
### 输出
- `output/` 中的修剪后 FASTQ 文件
- `output/trimmomatic.log`

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
- 依赖的 jar 和接头文件保留在 `dependencies/`。
- 适合作为整个转录组流程的第一步。
