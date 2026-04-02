# 8_2_coalescent/2_0_astral

## 功能概述
用于调用 ASTRAL 根据过滤后的 gene trees 推断 species tree。

## 正式入口
- `run_astral.py`

## 输入与输出
### 输入
- `input/` 中的 gene tree 文件
### 输出
- `output/` 中的 ASTRAL 物种树

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
- 支持可执行文件模式和 jar 模式。
