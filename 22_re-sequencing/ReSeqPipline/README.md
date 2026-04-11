# 22_re-sequencing/ReSeqPipline

## 功能概述
当前推荐用分步骤脚本运行流程，不再依赖单一主控脚本。

每个步骤都是一个独立的 Python 脚本，文件名前的数字就是流程顺序。每个脚本顶部都有“参数配置区”，直接改脚本里的变量即可，不需要命令行传参。

## 推荐入口
- `01_index_reference.py`
- `02_bwa_mapping.py`
- `03_gatk_call.py`
- `04_vcf_to_fasta.py`
- `05_check_non_atcg.py`
- `06_combine_fasta.py`

## 输入与输出
### 输入
- `input/reference/ref.fasta`
- `input/fastq/`

### 输出
- `output/bam_output/`
- `output/vcf_output/`
- `output/consensus_fasta_output/`
- `output/quality_report.txt`
- `output/combined_sequences.fasta`

## 使用方式
1. 先修改对应步骤脚本顶部的“参数配置区”。
2. 把输入文件放到 `input/` 下。
3. 按顺序运行需要的步骤脚本。
4. 调试时可以只跑某一个步骤。

## 说明
- `01` 到 `06` 每个脚本都内置了自己需要的环境变量和路径检查逻辑，可以单独运行。
- 当前流程脚本已经全部扁平化，主流程不再依赖 `dependencies/scripts/`。
