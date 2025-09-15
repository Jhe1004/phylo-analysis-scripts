# Trinity 转录组分析脚本工具包

本目录包含两个用于处理 Trinity 转录组组装结果和执行 Trinity 组装的 Python 脚本。以下是对这两个脚本的详细说明，包括其功能、使用方法和所需依赖。


---

## 1. `trinity_pe.py`

### 功能描述

该脚本用于自动化执行 Trinity 转录组组装流程，专门处理成对末端（Paired-End, PE）的测序数据。它会自动查找当前目录下符合命名规范的成对 `.fq` 或 `.fastq` 文件，并为每个样本生成并执行相应的 Trinity 命令。

- **输入**: 当前工作目录下，命名格式为 `<sample_name>_1.<ext>` 和 `<sample_name>_2.<ext>` 的成对测序文件（`.fq` 或 `.fastq`）。
- **输出**: 为每个样本在当前目录下创建一个名为 `<sample_name>_trinity.fasta` 的输出文件。
- **处理逻辑**:
  - 脚本首先扫描当前目录，查找所有符合 `_1` 和 `_2` 命名规则的成对文件。
  - 对于每对找到的文件，它会构建一个使用 Singularity 容器执行 Trinity 的命令。
  - 命令参数包括 (一般不需要用户自己设定)：
    - `--seqType fq`: 指定输入序列类型为 FASTQ。
    - `--left` 和 `--right`: 指定左端（_1）和右端（_2）的输入文件路径。
    - `--max_memory 40G`: 设置最大内存使用量为 40GB。
    - `--CPU 40`: 设置使用的 CPU 核心数为 40。
    - `--output`: 指定输出目录。
    - `--full_cleanup`: 在组装完成后进行清理。
  - 脚本会打印出为每个样本构建的完整命令，供用户检查。

### 使用方法

#### 安装前准备

1. **Python 环境**: 确保您的系统已安装 Python 3。
2. **Singularity**: 该脚本使用 Singularity 容器来运行 Trinity，因此需要预先安装 Singularity。
   - **安装 Singularity**:
     - **Linux**: 请参考 [Singularity 官方安装指南](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps)。通常需要管理员权限。
     - **Windows/macOS**: Singularity 主要为 Linux 设计。在 Windows 上，您可以使用 WSL2 (Windows Subsystem for Linux) 安装 Linux 发行版，然后在其中安装 Singularity。在 macOS 上，可以考虑使用虚拟机或容器服务。
3. **Trinity Singularity 镜像**: 脚本假定当前目录下存在名为 `trinityrnaseq.v2.15.1.simg` 的 Trinity Singularity 镜像文件。
   - 如果您没有此文件，需要从 Trinity 官方或通过 `singularity pull` 命令获取。例如：
     ```bash
     singularity pull docker://trinityrnaseq/trinityrnaseq:2.15.1
     ```
     这会生成 `trinityrnaseq_2.15.1.sif` 文件，您需要将其重命名为 `trinityrnaseq.v2.15.1.simg` 或修改脚本中的镜像名称。

#### 运行脚本

1. 确保所有成对的测序文件（如 `sampleA_1.fq`, `sampleA_2.fq`）都位于脚本所在的当前目录下。
2. 确保 `trinityrnaseq.v2.15.1.simg` 镜像文件也在当前目录下。
3. 打开终端或命令提示符（如果是 WSL2，请在 WSL2 环境中操作）。
4. 导航到包含此脚本、测序文件和 Singularity 镜像的目录。
5. 运行以下命令：
   ```bash
   python trinity_pe.py
   ```
6. 脚本会扫描目录，找到所有成对文件，并为每个样本打印出将要执行的 Trinity 命令。


### 配置说明

- **文件后缀名**: 脚本顶部的 `VALID_EXTENSIONS` 变量定义了要查找的文件后缀名列表。默认为 `['.fq', '.fastq']`。如果您使用的文件后缀名不同（如 `.fastq.gz`），请在此处修改。

### 注意事项

- 脚本假定所有输入文件都在当前工作目录下。
- 输出文件 `<sample_name>_trinity.fasta` 会直接在当前目录下创建。
- 脚本中设置的内存 (`40G`) 和 CPU (`40`) 参数是示例值，请根据您的实际硬件资源进行调整。
- 使用 Singularity 容器可以确保 Trinity 运行环境的一致性，避免依赖问题。
- 由于 Trinity 组装是一个计算密集型任务，执行时间可能较长，请耐心等待。

## 1. `get_longest_isoform_seq_per_trinity_gene.py`

### 功能描述

该脚本用于从 Trinity 组装结果中提取每个基因（gene）的最长转录本（isoform）序列。Trinity 通常会为每个基因生成多个转录本，这个脚本可以帮助用户筛选出每个基因的代表性序列（最长的那个），以便后续分析。

- **输入**: 包含多个 `.fasta` 文件的目录，每个文件包含 Trinity 组装的转录本序列。
- **输出**: 对于每个输入的 `.fasta` 文件，生成一个对应的 `.fas` 文件，其中包含该文件中每个基因的最长转录本序列。
- **处理逻辑**:
  - 脚本会遍历指定目录下的所有 `.fasta` 文件。
  - 对于每个文件，它会解析 Trinity 的序列 ID（如 `TRINITY_DN1000_c0_g1_i1`），从中提取基因 ID（如 `TRINITY_DN1000_c0_g1`）。
  - 对于每个基因 ID，它会比较所有相关的转录本长度，保留最长的一个。
  - 最终结果会按照转录本长度降序排列，并保存为 `.fas` 文件。

### 使用方法

#### 安装前准备

1. **Python 环境**: 确保您的系统已安装 Python 3。
2. **Biopython 库**: 该脚本依赖于 Biopython 库来处理 FASTA 文件。
   - 安装方法: 打开终端或命令提示符，运行以下命令：
     ```bash
     pip install biopython
     ```
     或者，如果您有多个 Python 版本，可能需要使用 `pip3`：
     ```bash
     pip3 install biopython
     ```

#### 运行脚本

1. 打开终端或命令提示符。
2. 导航到包含此脚本和 `.fasta` 文件的目录。
3. 运行以下命令：
   ```bash
   python get_longest_isoform_seq_per_trinity_gene.py <input_directory>
   ```
   - `<input_directory>`: 替换为包含 `.fasta` 文件的目录路径。
   - 例如，如果 `.fasta` 文件在当前目录下的 `fasta_files` 文件夹中，命令为：
     ```bash
   python get_longest_isoform_seq_per_trinity_gene.py fasta_files
   ```
4. 脚本执行后，会在 `<input_directory>` 中为每个 `.fasta` 文件生成一个同名但扩展名为 `.fas` 的输出文件。

### 注意事项

- 脚本会在开始时打印一条警告信息：“NOTE - longest transcript isn't always the best transcript! ... consider filtering based on relative expression support ...”。这提醒用户，最长的转录本不一定是最优的，根据表达量等信息进行过滤可能是更好的选择。
- 如果某个 `.fasta` 文件中没有找到有效的转录本，脚本会跳过该文件并不生成对应的输出文件。
- 脚本假定 Trinity 的序列 ID 遵循特定的格式（以 `_g` 或 `comp_c` 结尾后跟 `_` 和转录本信息），如果格式不同，可能需要修改正则表达式 `gene_id_pattern`。

