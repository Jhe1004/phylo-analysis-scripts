# 22_re-sequencing/ReSeqPipline

## 功能概述
这是重测序主流程总调度脚本。脚本会顺序调用索引参考、比对、变异检测、共识序列生成以及 FASTA 合并等依赖脚本，统一管理整个 resequencing 工作流。

## 正式入口
- `pipeline_master.py`

## 输入与输出
### 输入
- `input/` 中的参考基因组和原始测序数据
### 输出
- `output/` 中的各步结果、VCF、共识序列和汇总矩阵

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
- 底层依赖脚本保存在 `dependencies/scripts/`。
