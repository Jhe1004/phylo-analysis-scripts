# Trimmomatic 批量质控与文件整理脚本

## 1. 工具简介
这是一个用于批量处理二代测序（NGS）双端数据（Paired-End Data）的自动化 Python 脚本。

它的主要功能是：
1.  **自动识别**：在指定目录下自动寻找成对的 FASTQ 测序文件（如 `R1` 和 `R2`）。
2.  **批量质控**：调用 **Trimmomatic** 工具对每个样本进行去接头（Adapter）、去低质量碱基等操作。
3.  **自动清洗**：程序运行结束后，会自动**删除**不需要的非配对（Unpaired）中间文件，只保留高质量的配对文件。
4.  **规范命名**：将输出文件重命名为简洁的格式（如 `sample_1.fq`），方便后续分析软件调用。

---

## 2. 环境依赖

运行此脚本需要满足以下环境要求：

* **操作系统**：Linux, macOS, 或 Windows (需配置好 Java 环境)。
* **Python 版本**：**Python 3.6** 及以上版本。
    * 本脚本仅使用 Python 标准库（`os`, `sys`, `subprocess`），**无需安装** `pandas` 或 `biopython` 等第三方库。
* **外部软件**：
    * **Java Runtime Environment (JRE)**: 用于运行 Trimmomatic。
    * **Trimmomatic**: 必须拥有 `trimmomatic-0.40.jar` (或类似版本) 文件。

---

## 3. 输入/输出数据格式

### 输入数据 (Input)
脚本默认在**当前目录**下查找文件。文件必须是**双端测序数据**，且文件名后缀需要统一。

**示例文件结构**：
RawData/
├── SampleA_1.fq.gz   (正向序列)
├── SampleA_2.fq.gz   (反向序列)
├── SampleB_1.fq.gz
├── SampleB_2.fq.gz
└── trimmomatic_script.py (本脚本)

*在配置中，我们将后缀分别设置为 `1.fq.gz` 和 `2.fq.gz`。*

### 输出数据 (Output)
程序运行完毕后，目录中将生成质控后的文件。

**输出示例**：
RawData/
├── SampleA_1.fq      (质控后的 R1，已重命名)
├── SampleA_2.fq      (质控后的 R2，已重命名)
├── SampleB_1.fq
└── SampleB_2.fq

*注：生成的 `_unpaired.fq` 文件会被脚本自动删除。*

---

## 4. 使用方法

### 第一步：修改配置 (Configuration)
**无需在命令行输入任何参数**。请使用文本编辑器（如 VS Code, Notepad++, Sublime Text）打开脚本文件 `trimmomatic_script.py`。

在脚本顶部的 `### CONFIGURATION (配置区域) ###` 修改以下变量：

1.  `INPUT_DIRECTORY`: 数据的存放路径 (默认是当前路径)。
2.  `TRIMMOMATIC_JAR_PATH`: 你的 `trimmomatic-0.40.jar` 文件的实际路径。
3.  `INPUT_FORWARD_SUFFIX` / `INPUT_REVERSE_SUFFIX`: 你的原始文件后缀 (例如 `_R1.fastq.gz`)。
4.  `TRIM_STEPS`: 根据你的实验需求调整质控参数 (如 `SLIDINGWINDOW`, `MINLEN`)。

### 第二步：运行脚本
打开终端 (Terminal) 或命令行，进入脚本所在目录，输入以下命令：

python trimmomatic_script.py

脚本将自动开始处理，并在屏幕上打印进度日志。如果遇到错误（如找不到 Java），脚本会报错并提示原因。