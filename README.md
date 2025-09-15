# 系统发育基因组学分析脚本

本仓库包含了一套全面的系统发育基因组学分析脚本，涵盖了从原始序列数据处理到进化树构建和进化分析的完整流程。

## 目录
1. [概述](#概述)
2. [工作流程](#工作流程)
3. [目录结构和使用方法](#目录结构和使用方法)
   - [0_trimmomatic](#0_trimmomatic)
   - [1_trinity](#1_trinity)
   - [2_cdhit](#2_cdhit)
   - [3_TransDecoder-master](#3_transdecoder-master)
   - [4_1_proteinortho](#4_1_proteinortho)
   - [4_2_get_seq_from_proteinortho](#4_2_get_seq_from_proteinortho)
   - [5_1_del_chloro_mito_from_fasta](#5_1_del_chloro_mito_from_fasta)
   - [5_2_Treeshrink](#5_2_treeshrink)
   - [5_4_del_indel](#5_4_del_indel)
   - [6_1_concatenation](#6_1_concatenation)
   - [6_2_coalescent](#6_2_coalescent)
   - [7_daing](#7_daing)
   - [8_genetree_speciestree_conflict](#8_genetree_speciestree_conflict)
   - [11_Hybridization_detection](#11_hybridization_detection)
   - [12_WGD_detection](#12_wgd_detection)
   - [20_Chloroplast](#20_chloroplast)
   - [21_annotation](#21_annotation)
   - [22_re-sequencing](#22_re-sequencing)
   - [Other_tools](#other_tools)
4. [依赖要求](#依赖要求)
5. [使用指南](#使用指南)

## 概述

该分析流程专为全面的系统发育基因组学分析而设计，主要针对转录组数据。它包括数据预处理、组装、直系同源推断、序列比对、进化树构建以及杂交检测和全基因组复制分析等专门的进化分析工具。

## 工作流程

典型的工作流程包括以下步骤：
1. 原始测序数据的质量控制和修剪 (0_trimmomatic)
2. 转录组组装 (1_trinity)
3. 序列聚类和去除冗余 (2_cdhit)
4. 蛋白质预测 (3_TransDecoder-master)
5. 直系同源分析 (4_1_proteinortho, 4_2_get_seq_from_proteinortho)
6. 去除细胞器序列 (5_1_del_chloro_mito_from_fasta)
7. 系统发育树构建和优化 (5_2_Treeshrink, 5_4_del_indel, 6_1_concatenation)
8. 溯祖分析 (6_2_coalescent)
9. 专门的进化分析 (11_Hybridization_detection, 12_WGD_detection)

## 目录结构和使用方法

### 0_trimmomatic

使用Trimmomatic对原始双端测序数据进行质量控制和修剪的脚本。

- `trimmomatic_pe.py`: 处理双端FASTQ文件的主脚本。
  - 根据后缀自动检测配对文件
  - 运行Trimmomatic进行质量修剪
  - 清理未配对的reads并重命名输出文件
  - 在运行前配置后缀和Trimmomatic参数
  - 使用方法：将配对的FASTQ文件放在目录中，然后运行 `python trimmomatic_pe.py`

### 1_trinity

使用Trinity进行转录组组装的脚本。

- `trinity_pe.py`: 在双端reads上运行Trinity组装器。
  - 使用Singularity容器执行Trinity
  - 处理目录中的所有配对FASTQ文件
  - 直接在脚本中配置Trinity参数
  - 使用方法：将配对的FASTQ文件放在目录中，然后运行 `python trinity_pe.py`

### 2_cdhit

使用CD-HIT进行序列聚类的脚本和工具。

- 多个Perl脚本用于各种CD-HIT相关任务：
  - `cdhit.py`: 主CD-HIT执行脚本
  - 各种`clstr_*`脚本用于聚类分析和操作
  - `cd-hit-para.pl`: CD-HIT的并行执行
  - 使用方法：根据特定的聚类需求运行相应的脚本

### 3_TransDecoder-master

使用TransDecoder从转录组中预测蛋白质的脚本。

- `transdecoder.py`: 运行TransDecoder的主脚本。
  - 处理目录中的所有FASTA文件
  - 运行TransDecoder.LongOrfs和TransDecoder.Predict
  - 使用多进程进行并行执行
  - 使用方法：将转录组FASTA文件放在目录中，然后运行 `python transdecoder.py`

### 4_1_proteinortho

使用ProteinOrtho进行直系同源分析的脚本。

- `aln.py`: 直系同源分析中的序列比对脚本。

### 4_2_get_seq_from_proteinortho

根据ProteinOrtho结果提取序列的脚本。

- `get_seq_from_proteinortho.py`: 序列提取的主脚本。
- `heatmap_figure.py`: 生成直系同源结果的热图可视化。
- `violin_figure.py`: 生成直系同源结果的小提琴图可视化。

### 5_1_del_chloro_mito_from_fasta

从FASTA文件中删除叶绿体和线粒体序列的脚本。

- `del_chloro_mito_from_fasta.py`: 删除细胞器序列的主脚本。

### 5_2_Treeshrink

使用TreeShrink进行系统发育树优化的脚本。

- `compare_contents_two_folders.py`: 比较树文件和比对文件的一致性。
- 工作流程目录包含TreeShrink流程的每个步骤的脚本：
  - `0_mafft`: 序列比对
  - `1_del_indel`: 删除indel
  - `2_raxml`: 使用RAxML构建进化树
  - `3_combine_tree_with_name`: 合并树文件并添加名称
  - `4_Run_TreeShrink`: 运行TreeShrink优化
  - `5_filter_alignment_use_treeshrink_output`: 使用TreeShrink输出过滤比对

### 5_4_del_indel

从序列比对中删除indel的脚本。

- `delmissingsite.py`: 删除indel的主脚本。

### 6_1_concatenation

序列拼接和超矩阵构建的脚本。

- 工作流程目录包含拼接和建树的每个步骤的脚本：
  - `0_concat`: 序列拼接
  - `1_del_indel`: 删除indel
  - `2_codon`: 密码子分析
  - `20_score_tree`: 树评分
  - `3_raxml`: 使用RAxML构建进化树

### 6_2_coalescent

基于溯祖的物种树估计脚本。

- `2_0_astral/run_astral.py`: 运行ASTRAL溯祖分析的脚本。
- 工作流程目录包含溯祖分析的每个步骤的脚本：
  - `0_raxml`: 使用RAxML构建基因树
  - `1_filter`: 过滤基因树
  - `2_0_astral`: ASTRAL分析
  - `2_1_mp-est`: MP-EST分析
  - `3_coalescent_simulations`: 溯祖模拟
  - `4_SCOG_statistics`: SCOG统计

### 7_daing

密度树可视化的脚本。

- `densitree`: 包含生成密度树可视化的脚本。

### 8_genetree_speciestree_conflict

分析基因树和物种树之间冲突的脚本。

- 包含QuIBL等基因树-物种树协调工具。

### 11_Hybridization_detection

使用多种方法检测杂交事件的脚本。

- `Hyde`: 使用系统发育网络进行杂交检测的脚本
  - `0_simulate_sequences`: 序列模拟
  - `1_run_hyde`: 运行HyDe分析
  - `2_visual_hyde`: HyDe结果可视化
- `Phylonet`: 使用PhyloNet构建系统发育网络的脚本
  - `0_将想要检测的树提取出来`: 提取待检测的树
  - `1_Get_phylonet_input`: 生成PhyloNet输入
  - `2_Run_phylonet`: 运行PhyloNet
  - `3_phylonet_result_to_dendroscope`: 将结果转换为Dendroscope格式
- `Phylonetworks`: 系统发育网络分析的脚本
  - `1_将想要检测的树提取出来`: 提取待检测的树
  - `2_Astral计算先导树`: 使用ASTRAL计算引导树
  - `3_phylonetworks计算`: PhyloNetworks计算
  - `4_phylonetworks_bootstrap`: PhyloNetworks自举分析
  - `5_phylonetworks_figure`: PhyloNetworks结果可视化
- `QuIBL`: 使用贝叶斯线性模型量化基因渐渗的脚本
  - `QuIBL.py`: 主程序
  - `visual_quibl.py`: 结果可视化
  - `analysis`: 分析相关脚本
  - `cython_vers`: Cython版本
  - `Small_Test_Example`: 小测试示例

### 12_WGD_detection

全基因组复制检测的脚本。

- `absolute_dating`: 时间校准的系统发育分析脚本
  - `0_get_gene_family`: 获取基因家族
  - `1_mafft`: 序列比对
  - `2_raxml`: 使用RAxML构建基因树
  - `3_reroot_gene_family_trees`: 重新根植基因家族树
  - `4_generate_treePL_file`: 生成TreePL文件
  - `6_extract_time`: 提取时间信息
  - `7_result`: 结果分析
- `ks_based`: 基于Ks的WGD检测脚本
  - `0_selfblast`: 自比对
  - `1_calculate_KS_WGD`: 计算WGD的Ks值
  - `2_calculate_KS_two_species`: 计算两个物种间的Ks值
  - `new`: 新方法
- `phypart`: 基因树分区分析脚本
  - `0_get_gene_family`: 获取基因家族
  - `1_mafft`: 序列比对
  - `2_del_indel`: 删除indel
  - `3_raxml`: 使用RAxML构建基因树
  - `4_reroot_gene_family_trees`: 重新根植基因家族树
  - `5_rename`: 重命名
  - `6_phypart`: PhyPart分析

### 20_Chloroplast

叶绿体基因组分析的脚本。

- `0_密码子偏好`: 密码子使用偏好分析
  - `1_extract_cds.py`: 提取CDS序列
  - `2_calculate_rscu.py`: 计算RSCU值
  - `3_plot_heatmap_with_tree.py`: 绘制带树的热图
- `1_ssr分析`: SSR(简单序列重复)分析
  - `1_gb_to_fasta.py`: GenBank格式转换为FASTA格式
  - `2_find_ssrs.py`: 查找SSRs
  - `3_summarize_and_plot_ssrs.py`: 总结和绘制SSRs
- `2_基因正选择分析`: 正选择分析
- `3_基因表格生成`: 基因表格生成
- `4_高变区`: 高变区识别
- `5_mvista`: 基于mvista的分析
- `6_IRscope`: 反向重复分析
- `chloroplast_genome_alignment-master`: 叶绿体基因组比对工具
- `PGA`: 叶绿体基因组组装工具

### 21_annotation

基因注释的脚本。

- `annotate_fasta_genes`: 基因注释工具
- `Genoma`: 基因组注释工具

### 22_re-sequencing

重测序数据分析的脚本。

- `GetGeneFromGFF`: 从GFF文件中提取基因序列
- `ReSeqPipline`: 完整的重测序分析流程

### Other_tools

其他生物信息学工具。

- `0_树相关`: 树相关的工具
  - `0_删掉树中的一些类群`: 从系统发育树中删除指定类群
  - `1_提取树中的所有物种名称`: 提取系统发育树中的所有物种名称
  - `2_重命名树中的所有样品名称`: 重命名系统发育树中的样品名称
  - `3_对树重新置根`: 对系统发育树重新置根
  - `4_输入物种名称随机生成进化树`: 根据输入的物种名称随机生成系统发育树
- `1_统计各种数字`: 统计工具
  - `0_fasta序列N50统计`: 统计FASTA序列的N50值
  - `1_统计get_proteinortho结果中各个物种的基因数目`: 统计ProteinOrtho结果中各物种的基因数
  - `3_alignment整体信息位点统计`: 统计序列比对的整体信息位点
  - `4_alignment各序列gap信息统计`: 统计序列比对中各序列的gap信息
- `2_画各种图`: 可视化工具
  - `0_频率分布直方图`: 绘制频率分布直方图
  - `1_小提琴图`: 绘制小提琴图
- `3_fasta&fastq相关`: FASTA/FASTQ处理工具
  - `0_重命名fasta中的所有样品名称`: 重命名FASTA文件中的样品名称
  - `1_bootstrap_fasta`: 对FASTA序列进行bootstrap抽样
  - `2_提取fasta矩阵中的部分位点`: 从FASTA矩阵中提取部分位点
  - `3_提取fasta矩阵中的部分序列`: 从FASTA矩阵中提取部分序列
  - `4_将一个大的fastq文件分割开`: 将大的FASTQ文件分割成小文件

## 依赖要求

- Python 3.x
- Java (用于Trimmomatic, ASTRAL等工具)
- Trinity (转录组组装器)
- CD-HIT (序列聚类)
- TransDecoder (蛋白质预测)
- RAxML (系统发育树构建)
- ASTRAL (溯祖分析)
- TreeShrink (树优化)
- PhyloNet (系统发育网络构建)
- PhyloNetworks (系统发育网络分析)
- 其他在特定脚本中引用的各种生物信息学工具

## 使用指南

1. 按照工作流程顺序操作以获得最佳结果
2. 根据您的系统和数据在每个脚本中配置路径和参数
3. 确保所有依赖项都已正确安装并可访问
4. 查看各个脚本的注释了解具体使用说明
5. 一些脚本使用容器(Singularity/Docker) - 确保这些容器已正确配置
6. 对于复杂的工作流程(如杂交检测和WGD检测)，按照子目录中的数字顺序执行脚本