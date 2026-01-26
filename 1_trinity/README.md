# 1. Trinity 组装 & 冗余提取流水线 🧬

本目录包含了处理转录组原始数据（Raw Reads）到生成非冗余组装序列（Non-redundant Trinity Assembly）的完整自动化方案。

## 🌟 核心功能
1. **自动化批量组装 (`01_batch_trinity.py`)**: 
   - 自动识别双端测序文件。
   - 调用 Singularity 容器运行 Trinity，无需在宿主机安装复杂的生物信息依赖。
   - 自动清理中间文件，极大节省磁盘空间。
2. **最长转录本提取 (`02_extract_longest_isoform.py`)**: 
   - 针对 Trinity 产生的大量 Isoforms，采用贪婪算法保留每个基因最长的序列。
   - 显著降低后续系统发育分析的计算压力。

---

## 🚀 快速开始

### 1. 下载必要的镜像文件 (Singularity Image)
由于镜像文件（约 3GB）超过了 GitHub 的限制，请从 Hugging Face 下载：
- **下载地址**: [Hugging Face - phylo-scripts-images](https://huggingface.co/datasets/Jhe1004/phylo-scripts-images/resolve/main/trinityrnaseq.v2.15.2.simg)
- **存放位置**: 下载后请将 `trinityrnaseq.v2.15.2.simg` 放在本目录下。

### 2. 环境要求
- **Python**: 3.6+ (需安装 `biopython` 库)
- **容器软件**: 系统需安装 [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) 或 Apptainer。

```bash
pip install biopython
```

---

## 📖 详细步骤

### 第一步：批量组装 (01_batch_trinity.py)
1. **配置**: 修改脚本顶部的 `CONFIGURATION` 区域（输入目录、线程数、内存上限）。
2. **测试**: 默认 `DRY_RUN = True`，运行后仅打印命令供检查。
3. **执行**:
   ```bash
   python 01_batch_trinity.py
   ```

### 第二步：提取最长转录本 (02_extract_longest_isoform.py)
该脚本会扫描已完成的组装目录，提取 `Trinity.fasta` 中的最长转录本。
1. **执行**:
   ```bash
   python 02_extract_longest_isoform.py
   ```
2. **输出**: 每个样本会生成一个 `.longest.fas` 文件，用于后续分析。

---

## 📂 推荐目录结构
```text
1_trinity/
├── 01_batch_trinity.py              # 组装脚本
├── 02_extract_longest_isoform.py    # 提取脚本
├── trinityrnaseq.v2.15.2.simg       # 从 Hugging Face 下载的镜像
├── raw_data/                        # [输入] 原始测序数据 (.fq.gz)
└── assembly_results/                # [输出] 组装结果
```

---

## ⚠️ 注意事项
- **内存建议**: Trinity 组装非常耗费内存，建议每个样本分配 50GB 以上内存。
- **清除率**: 提取最长转录本后，序列总数通常会减少 40%-60%，属于正常现象（去除了冗余的异构体）。
