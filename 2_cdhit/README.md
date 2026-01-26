# 2. CD-HIT 并行去冗余工具 ✂️

本目录提供了一个高效的并行处理脚本，用于利用 **CD-HIT-EST** 工具对组装好的转录组序列进行去冗余处理。

## 🌟 核心功能
- **Singularity 容器化**: 无需在本地安装 CD-HIT 环境，直接调用镜像运行。
- **全自动并行**: 自动检测 CPU 核心数，同时对目录下所有 `.fas` 文件进行批量处理。
- **高兼容性**: 专为大规模转录组分析设计，能够处理海量序列数据。

---

## 🚀 快速开始

### 1. 下载必要的镜像文件 (Singularity Image)
由于镜像文件较大，请从 Hugging Face 下载：
- **下载地址**: [Hugging Face - phylo-scripts-images](https://huggingface.co/datasets/Jhe1004/phylo-scripts-images/resolve/main/cdhit.sif)
- **存放位置**: 下载后请将 `cdhit.sif` 放在本目录下。

### 2. 环境要求
- **Python**: 3.x
- **容器软件**: 系统需安装 [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) 或 Apptainer。

---

## 📖 使用指南

### 运行程序
1. 将需要处理的 `.fas` 文件（如上一步得到的最长转录本）放入本目录，或在 `cdhit.py` 脚本中指定路径。
2. 执行脚本：
   ```bash
   python cdhit.py
   ```

### 参数说明
脚本默认使用以下 CD-HIT 参数：
- `identity = 0.98` (序列相似度阈值)
- `coverage = 0.9` (序列覆盖度阈值)
- 更多详细参数可进入 `cdhit.py` 内部简单修改。

---

## 📂 输出文件说明
处理完成后，每个输入文件 `sample.fas` 会对应产生：
1. `sampleta`: (默认输出名后缀) 聚类去冗余后的非冗余序列文件。
2. `sampleta.clstr`: 聚类详细信息记录，包含哪些序列归为了同一类。
3. `processing.log`: 记录运行过程中的并行详细日志。

---

## 🛠 故障排除
- **镜像丢失**: 如果报错 `Image not found`，请确认 `cdhit.sif` 已经放置在当前脚本目录下。
- **权限问题**: 如果 Singularity 运行失败，请尝试赋予镜像执行权限：`chmod +x cdhit.sif`。