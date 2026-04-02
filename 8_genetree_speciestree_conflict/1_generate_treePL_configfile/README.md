# 8_genetree_speciestree_conflict/1_generate_treePL_configfile

## 功能概述
用于根据树和定年约束批量生成 treePL 配置文件，为大量 gene tree 的定年分析做准备。

## 正式入口
- `generate_treePL_conf.py`

## 输入与输出
### 输入
- `input/` 中的树文件和校准信息
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
