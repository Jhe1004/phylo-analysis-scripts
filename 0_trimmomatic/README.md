# Trimmomatic 脚本

此目录包含用于自动化处理 Trimmomatic 预处理测序数据常见任务的 Python 脚本。

## 脚本

### 1. `gunzip.py`

此脚本自动化解压缩当前目录中的 `.gz` 文件。

#### 依赖

- `gunzip` 工具（大多数类 Unix 系统的一部分）
  - 注意：在 Windows 上，您可能需要安装额外的工具，如 Git Bash 或 WSL 才能使用 `gunzip`。

#### 配置

在运行脚本之前，您可能需要修改脚本中的以下参数：

- `TARGET_EXTENSION`：要查找的文件扩展名（默认为 `.gz`）。
- `TARGET_DIRECTORY`：搜索文件的目录（默认为当前目录 `.`）。

#### 使用方法

1. 将脚本放置在包含要解压缩的 `.gz` 文件的目录中。
2. 运行脚本：
   ```bash
   python gunzip.py
   ```

#### 输入文件

- 目标目录中具有指定 `TARGET_EXTENSION` 的任何文件。

#### 输出文件

- 脚本将就地解压缩 `.gz` 文件，删除原始压缩文件并创建未压缩的版本。

### 2. `trimmomatic_pe.py`

此脚本自动化在双端 FASTQ 文件上运行 Trimmomatic 的过程。

#### 依赖

- Java（用于运行 Trimmomatic）
- Trimmomatic JAR 文件（`trimmomatic-0.40.jar` 应该在同一个目录中）

#### 配置

在运行脚本之前，您可能需要修改脚本中的以下参数：

- `forward_suffix`：正向读取文件的后缀（默认为 `1.fq.gz`）。
- `reverse_suffix`：反向读取文件的后缀（默认为 `2.fq.gz`）。
- `trimmomatic_jar`：Trimmomatic JAR 文件的路径（默认为 `trimmomatic-0.40.jar`）。
- `threads`：要使用的线程数（默认为 `30`）。
- `phred`：质量编码（`-phred33` 或 `-phred64`，默认为 `-phred33`）。
- `trim_options`：修剪参数（默认为 `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36`）。

#### 使用方法

1. 将您的双端 FASTQ 文件放置在目录中。
2. 确保 `trimmomatic-0.40.jar` 在同一个目录中。
3. 如果需要，修改脚本中的配置参数。
4. 运行脚本：
   ```bash
   python trimmomatic_pe.py
   ```

#### 输入文件

- 具有指定 `forward_suffix` 和 `reverse_suffix` 的双端 FASTQ 文件。
- 适配器文件（如果在修剪中使用）应该在 `adapters/` 目录中。

#### 输出文件

- `<sample>_1_paired.fq` 和 `<sample>_2_paired.fq`：通过质量过滤的配对读取。
- `<sample>_1_unpaired.fq` 和 `<sample>_2_unpaired.fq`：通过质量过滤的未配对读取。
- 处理后，未配对的文件将被删除，配对的文件将重命名为 `<sample>_1.fq` 和 `<sample>_2.fq`。