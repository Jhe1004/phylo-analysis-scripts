# Transcriptome Phylogenomics Pipeline

本仓库包含了一套统一的系统发育基因组学分析脚本。所有脚本均已标准化，通过命令行参数运行。

## 目录结构
- `lib/`: 共享工具库
- `0_trimmomatic/`: 质控与修剪
- `1_trinity/`: 转录组组装
- `2_cdhit/`: 序列去冗余
- (其他目录保持不变)

## 通用用法
所有脚本现在都使用标准命令行参数。您不再需要编辑脚本内部变量。
运行 `python script.py --help` 查看详细帮助。

### 0. 质控 (Trimmomatic)
```bash
cd 0_trimmomatic
python trimmomatic_pe.py \
    --input /path/to/raw_data \
    --jar /path/to/trimmomatic.jar \
    --threads 30
```
- **自动检测**: 默认查找 `_1.fq.gz` 和 `_2.fq.gz` 结尾的配对文件。
- **参数**: 可通过 `--trim-params` 自定义修剪参数。

### 1. 组装 (Trinity)
```bash
cd 1_trinity
python 01_batch_trinity.py \
    --input /path/to/clean_data \
    --threads 40 \
    --memory 50G
```
- **自动检测**: 支持 `.fq.gz`, `.fastq`, `.fasta` 等多种格式。

### 2. 去冗余 (CD-HIT)
```bash
cd 2_cdhit
python cdhit.py \
    --input /path/to/assemblies \
    --threshold 0.95 \
    --threads 10
```
- **输入**: 自动处理目录下所有 `.fas` 文件。

## 依赖
- Python 3+
- Trimmomatic (Java)
- Trinity
- CD-HIT