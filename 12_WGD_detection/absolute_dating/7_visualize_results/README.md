# absolute_dating Step 7

## 功能概述
把提取出的年龄结果绘制成密度图、热图和综合统计图，用于识别潜在 WGD 峰。

## 正式入口
- `run_result_plots.py`

## 输入与输出
### 输入
- `input/` 中的 sample ages、species tree、分组信息
### 输出
- `output/` 中的图像和统计表

## 使用方式

1. 打开正式脚本，直接修改脚本最前面的配置区。
2. 把输入文件放入当前目录的 `input/`。
3. 在当前目录运行对应脚本。
4. 结果会统一写入 `output/`。

本项目已经统一约定：

- 不通过 CLI 交互式传参
- 路径尽量写相对路径
- 外部软件优先从 `dependencies/bin/`、指定 conda 环境、系统 `PATH` 查找
