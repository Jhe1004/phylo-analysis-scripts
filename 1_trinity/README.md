# Trinity 组装自动化与最长转录本提取流程 (Auto-Trinity Pipeline)

# 1. 项目背景与工具简介

在植物转录组学（Transcriptomics）研究中，使用 Trinity 进行 De novo 组装是标准流程之一。然而，面对几十甚至上百个样本时，手动编写脚本容易出错且效率低下。此外，Trinity 组装结果通常包含大量的冗余异构体（Isoforms），直接用于后续的系统发育分析（Phylogenetics）会增加计算负担且引入噪音。

本项目提供了一套标准化的 Python (v3.6+) 自动化流水线，包含两个核心脚本：

1. **01_batch_trinity.py**: 基于 Singularity 容器的 Trinity 批量组装工具。

2. **02_extract_longest_isoform.py**: 基于 Biopython 的基因代表序列提取工具（清洗冗余 Isoform）。

---

# 2. 脚本工作原理 (Algorithm & Logic)

为了确保分析结果的可信度，了解脚本背后的处理逻辑至关重要。

## 脚本 1: 批量组装 (01_batch_trinity.py)

该脚本采用 "扫描-配对-执行" 的逻辑：

1. **文件扫描**: 遍历指定目录，查找所有符合后缀（如 .fq.gz）的文件。

2. **智能配对**:
        脚本依据文件名中的 _1 和 _2 标识自动识别双端测序数据（Paired-End Reads）。

3. 只有当 Sample_1.fq 和 Sample_2.fq 同时存在时，才会被视为有效样本。

4. **命令构建**:
        自动为每个样本生成独立的输出目录（防止文件冲突）。

5. 基于用户配置（CPU/内存），动态构建调用 Singularity 容器的 Trinity 命令。

6. 命令包含 --full_cleanup 参数，确保运行结束后自动删除庞大的中间临时文件，节省硬盘空间。

## 脚本 2: 提取最长转录本 (02_extract_longest_isoform.py)

该脚本采用 "流式读取-字典更新" 的贪婪算法（Greedy Algorithm）：

1. **流式解析**: 使用 Biopython.SeqIO 逐行读取巨大的 FASTA 文件，内存占用低。

2. **ID 解析**: 利用正则表达式（RegEx）从复杂的 Trinity ID（如 TRINITY_DN100_c0_g1_i1）中提取核心 Gene ID（TRINITY_DN100_c0_g1）。

3. **长度比较**:
        维护一个字典 Dict[GeneID, SequenceRecord]。

4. 对于每个读入的序列，检查其 Gene ID 是否已存在于字典中。

5. **逻辑**: 如果是新基因，直接存入；如果基因已存在，比较当前序列与字典中序列的长度，保留较长的那一条。

6. **排序输出**: 处理完所有序列后，将结果按序列长度从长到短排序输出，便于后续观察统计。

---

# 3. 环境依赖 (Dependencies)

## 基础环境

- **操作系统**: Linux (CentOS 7+, Ubuntu 18.04+)

- **Python**: 3.6 或更高版本 (必须支持 f-string 和 pathlib)

- **容器引擎**: Singularity (用于运行 Trinity，避免复杂的依赖配置)

## Python 库安装

脚本 2 依赖 Biopython 进行序列处理，请使用 pip 安装：

```bash

pip install biopython
```

## 外部文件

Trinity 镜像: 必须拥有 Trinity 的 Singularity 镜像文件（.simg 或 .sif）。

推荐版本: trinityrnaseq.v2.15.2

---

# 4. 输入与输出数据详解

## 目录结构示例

假设您的工作目录如下：

```plaintext

Project_Root/
├── 01_batch_trinity.py <-- 脚本 1
├── 02_extract_longest_isoform.py <-- 脚本 2
├── trinityrnaseq.v2.15.2.simg <-- 镜像文件
├── raw_data/ <-- [输入] 原始数据目录
│   ├── SampleA_1.fq.gz
│   ├── SampleA_2.fq.gz
│   ├── SampleB_1.fq.gz
│   └── SampleB_2.fq.gz
└── assembly_results/ <-- [输出] 结果将生成在这里
```

## 脚本 1 的数据流

- 输入: raw_data/ 中的双端压缩文件。

- 格式要求: 必须包含 _1 和 _2 结尾。

- 输出: 脚本会在当前目录下生成每个样本的文件夹。
        SampleA_trinity/Trinity.fasta

- SampleB_trinity/Trinity.fasta

## 脚本 2 的数据流

- 输入: 上一步生成的 Trinity.fasta 文件。

- ID 格式示例: >TRINITY_DN1000_c0_g1_i1 len=500 path=[...]

- 输出: 同名目录下生成的 .longest.fas 文件。
        文件内容: 每个 Gene 只保留一条序列。

- ID 格式保持不变，但文件体积通常减少 30%-50%。

---

# 5. 详细使用指南

注意: 本工具采用“配置驱动”模式。无需记忆复杂的命令行参数，直接修改脚本顶部的配置区即可。

## 步骤一：批量组装 (01_batch_trinity.py)

### 编辑配置:

用文本编辑器（如 vim 或 nano）打开脚本，修改 ### CONFIGURATION ### 区域：

```python

# CONFIGURATION ###
INPUT_DIRECTORY = "./raw_data"  # 你的测序文件路径
VALID_EXTENSIONS = ['.fq.gz', '.fastq.gz']
SINGULARITY_IMAGE_PATH = "./trinityrnaseq.v2.15.2.simg"  # 镜像路径
CPU_THREADS = 40  # 建议设为服务器最大核心数
MAX_MEMORY = "100G"  # Trinity 非常吃内存，建议 50G+
DRY_RUN = True  # 首次运行建议设为 True，检查命令是否正确
```

### 执行:

```bash

python 01_batch_trinity.py
```

如果 DRY_RUN = True，屏幕将打印出所有生成的 Trinity 命令。确认无误后，将 DRY_RUN 改为 False 并重新运行。

## 步骤二：提取最长转录本 (02_extract_longest_isoform.py)

### 准备数据:

你可以将所有组装好的 Trinity.fasta 移动到一个单独的文件夹（例如 assembly_results），或者直接在当前目录下运行。

### 编辑配置:

打开脚本修改配置：

```python

# CONFIGURATION ###
INPUT_DIRECTORY = "./assembly_results"  # 包含 Trinity.fasta 的文件夹
INPUT_FILE_PATTERN = "*.fasta"  # 匹配所有 fasta 文件
OUTPUT_EXTENSION = ".longest.fas"  # 输出后缀
```

### 执行:

```bash

python 02_extract_longest_isoform.py
```

脚本会输出处理日志：

```plaintext

正在读取: SampleA.fasta ...
-> 完成。输入序列: 45000, 提取基因: 23000 (过滤率 48.9%)
```

---

# 6. 常见问题 (FAQ)

### Q1: 脚本 1 报错 exec: "singularity": executable file not found

A: 您的服务器未安装 Singularity，或者没有将其添加到环境变量 $PATH 中。请联系管理员安装。

### Q2: 为什么提取后的序列数量变少了？

A: 这是预期的结果。Trinity 会为同一个基因组装出多个剪接变体（Isoform，例如 i1, i2, i3）。本脚本只保留最长的一条（通常代表最具代表性的编码区），因此序列总数会下降。


