# GetGeneFromGFF Step 0

## 功能概述
从原始 GFF 中提取 mRNA 或目标转录本相关注释。

## 正式入口
- `extract_mRNA_annotation.py`

## 输入与输出
### 输入
- `input/` 中的 GFF 文件
### 输出
- `output/` 中的精简注释文件

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
