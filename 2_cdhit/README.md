# CD-HIT 并行处理脚本

这是一个用于并行处理 FASTA 文件的 Python 脚本，使用 CD-HIT-EST 工具对转录组序列进行聚类分析。

## 功能描述

- **批量处理**: 自动扫描当前目录下所有 `.fas` 文件
- **并行处理**: 使用多进程并行处理多个文件，提高效率
- **自动检测**: 自动检测系统中安装的 CD-HIT-EST 可执行文件
- **日志记录**: 详细的日志记录，包括处理进度和错误信息

## 系统要求

- Python 3.x
- CD-HIT 软件包 (cd-hit-est 可执行文件)

## 安装 CD-HIT

### Ubuntu/Debian
```bash
sudo apt-get install cd-hit
```

### CentOS/RHEL
```bash
sudo yum install cd-hit
```

### 源码编译
```bash
# 从官网下载源码
wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz

# 解压并编译
tar -xzf cd-hit-v4.8.1-2019-0228.tar.gz
cd cd-hit-v4.8.1-2019-0228
make

# 将可执行文件添加到 PATH
sudo cp cd-hit-est /usr/local/bin/
```

## 使用方法

### 基本用法
```bash
python cdhit.py
```

### 指定 CD-HIT-EST 路径
如果 CD-HIT-EST 不在系统路径中，可以指定完整路径：
```bash
python cdhit.py /path/to/cd-hit-est
```

## 输入文件要求

- 文件格式: `.fas` 后缀的 FASTA 文件
- 文件位置: 脚本所在目录
- 序列类型: 转录组序列 (EST)

## 输出文件

- 每个输入文件会生成对应的输出文件，命名格式为 `原文件名ta`
- 处理日志保存在 `processing.log` 文件中

## 输出文件说明

CD-HIT-EST 会生成以下文件：
- `*.clstr`: 聚类结果文件，包含序列聚类信息
- 主输出文件: 去冗余后的序列文件

## 并行处理

脚本会自动检测 CPU 核心数，并使用所有可用核心并行处理文件。

## 日志信息

脚本会输出详细的处理信息，包括：
- 找到的 `.fas` 文件数量
- 并行进程数量
- 每个文件的处理状态
- 错误信息（如果有）

## 示例

假设当前目录有以下文件：
```
sample1.fas
sample2.fas
sample3.fas
```

运行脚本：
```bash
python cdhit.py
```

输出：
```
2024-01-01 12:00:00 [INFO] Found 3 .fas files to process.
2024-01-01 12:00:00 [INFO] Starting processing with 8 parallel processes...
2024-01-01 12:00:00 [INFO] Starting cd-hit-est for sample1.fas -> sample1.fasta
2024-01-01 12:00:00 [INFO] Starting cd-hit-est for sample2.fas -> sample2.fasta
2024-01-01 12:00:00 [INFO] Starting cd-hit-est for sample3.fas -> sample3.fasta
2024-01-01 12:00:05 [INFO] Successfully processed sample1.fas
2024-01-01 12:00:06 [INFO] Successfully processed sample2.fas
2024-01-01 12:00:07 [INFO] Successfully processed sample3.fas
2024-01-01 12:00:07 [INFO] All files have been processed.
```

## 故障排除

### CD-HIT-EST 未找到
如果出现 "cd-hit-est 未安装或不在系统路径中" 错误：
1. 确保 CD-HIT 已正确安装
2. 或者使用完整路径运行脚本：`python cdhit.py /path/to/cd-hit-est`

### 没有找到 .fas 文件
确保当前目录包含 `.fas` 后缀的 FASTA 文件。

### 权限问题
确保 CD-HIT-EST 可执行文件具有执行权限：
```bash
chmod +x /path/to/cd-hit-est
```

## 注意事项

- 脚本会处理当前目录下的所有 `.fas` 文件
- 输出文件会覆盖同名的现有文件
- 建议在处理前备份重要数据
- 处理大型文件时可能需要较长时间和较多内存