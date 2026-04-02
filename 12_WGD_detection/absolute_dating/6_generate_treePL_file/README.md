# absolute_dating Step 6

## 功能概述
批量生成 treePL 输入配置文件，把基因树、校准信息和输出路径统一组织好。

## 正式入口
- `generate_treePL_file.py`

## 输入与输出
### 输入
- `input/` 中的 gene trees 与校准信息
### 输出
- `output/` 中的 treePL 配置文件

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
