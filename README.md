# transcriptome_pyloneny

## 项目概述

这是一个围绕转录组系统发育分析而逐步扩展出来的脚本仓库，后续又加入了叶绿体分析、杂交检测、WGD 检测、注释、重测序以及一批常用辅助工具。  
这个根目录 `README.md` 的目标不是替代各子目录自己的说明，而是帮读者先建立一张“仓库地图”，快速知道：

- 每一类分析放在哪个文件夹
- 每个文件夹里的正式入口脚本是什么
- 多步骤流程的顺序大致如何
- 遇到 `README.md`、`dependencies/`、上游工具时应该怎么看

## 仓库整理原则

项目中的大部分流程目录都尽量遵循下面的约定：

- 正式运行脚本放在流程目录顶层或步骤子目录顶层
- 输入数据尽量放在 `input/`
- 输出结果尽量写到 `output/`
- 用户需要修改的参数尽量集中在脚本开头配置区
- 依赖脚本、模板文件、第三方上游程序一般放在 `dependencies/`
- 如果某个目录下有总流程，通常会同时提供一个同级 `README.md`

## 建议阅读顺序

1. 先看本文件，确认目标分析属于哪个大目录
2. 再进入对应目录阅读该目录自己的 `README.md`
3. 如果该目录下面还有分步骤子目录，再按编号顺序阅读子目录 `README.md`
4. 最后修改正式脚本顶部配置区并运行

## 根目录总览

### 转录组主流程目录

- `0_trimmomatic/`：原始双端测序数据质控与修剪
- `1_trinity/`：Trinity 组装与最长转录本提取
- `2_cdhit/`：CD-HIT 去冗余聚类
- `3_TransDecoder-master/`：ORF 预测、蛋白与 CDS 提取
- `4_1_proteinortho/`：同源基因簇识别
- `4_2_get_seq_from_proteinortho/`：从 proteinortho 结果提取序列
- `5_del_chloro_mito_from_fasta/`：去除叶绿体和线粒体相关序列
- `6_1_mafft/`：批量比对
- `6_2_del_indel/`：缺失位点和 indel 过滤
- `7_Treeshrink/`：建树并用 TreeShrink 清理异常长枝

### 系统发育树与定年相关目录

- `8_1_concatenation/`：拼接矩阵建树流程
- `8_2_coalescent/`：共祖模型建树流程
- `8_genetree_speciestree_conflict/`：基因树与物种树冲突分析
- `9_daing/`：treePL 单次定年

### 杂交与 WGD 检测

- `11_Hybridization_detection/`：杂交检测
- `12_WGD_detection/`：WGD 检测与定年相关分析

### 专题分析目录

- `20_Chloroplast/`：叶绿体专题分析
- `21_annotation/`：注释相关流程
- `22_re-sequencing/`：重测序流程和基因提取流程

### 辅助目录

- `Other_tools/`：零散常用工具脚本
- `lib/`：仓库内部通用 Python 工具函数
- `脚本整理规范.md`：脚本整理时采用的规范说明

## 详细目录地图

下面按目录列出正式脚本入口。默认不展开 `dependencies/` 中的上游程序和辅助脚本，因为那些通常不是推荐的直接入口。

### 0_trimmomatic

- 作用：原始双端 reads 质控、剪切接头、过滤低质量序列
- 正式脚本：
- `0_trimmomatic/trimmomatic_pe.py`
- 补充：
- `0_trimmomatic/README.md`
- `0_trimmomatic/dependencies/` 中保留了 `trimmomatic-0.40.jar` 和 adapter 文件

### 1_trinity

- 作用：批量 Trinity 组装，并提取最长 isoform
- 正式脚本：
- `1_trinity/trinity_pipeline.py`
- 补充：
- `1_trinity/README.md`
- `1_trinity/dependencies/scripts/` 中有流程辅助脚本

### 2_cdhit

- 作用：转录本或蛋白序列去冗余聚类
- 正式脚本：
- `2_cdhit/cdhit.py`
- 补充：
- `2_cdhit/README.md`

### 3_TransDecoder-master

- 作用：ORF 预测并输出蛋白、CDS 等结果
- 正式脚本：
- `3_TransDecoder-master/transdecoder.py`
- 补充：
- `3_TransDecoder-master/README.md`
- `3_TransDecoder-master/dependencies/upstream/` 中保留了上游 TransDecoder 工具文件

### 4_1_proteinortho

- 作用：跨物种同源基因簇识别
- 正式脚本：
- `4_1_proteinortho/run_proteinortho.py`
- 补充：
- `4_1_proteinortho/README.md`

### 4_2_get_seq_from_proteinortho

- 作用：从 proteinortho 结果中提取 pep / cds 序列
- 正式脚本：
- `4_2_get_seq_from_proteinortho/get_seq_from_proteinortho.py`
- 补充：
- `4_2_get_seq_from_proteinortho/README.md`

### 5_del_chloro_mito_from_fasta

- 作用：去除疑似叶绿体或线粒体来源序列
- 正式脚本：
- `5_del_chloro_mito_from_fasta/del_chloro_mito_from_fasta.py`
- 补充：
- `5_del_chloro_mito_from_fasta/README.md`

### 6_1_mafft

- 作用：批量多序列比对
- 正式脚本：
- `6_1_mafft/mafft.py`
- 补充：
- `6_1_mafft/README.md`

### 6_2_del_indel

- 作用：删除缺失比例高或含 indel 的位点
- 正式脚本：
- `6_2_del_indel/delmissingsite.py`
- 补充：
- `6_2_del_indel/README.md`

### 7_Treeshrink

- 作用：先建树，再运行 TreeShrink 清理异常长枝序列
- 正式脚本：
- `7_Treeshrink/run_pipeline.py`
- 补充：
- `7_Treeshrink/README.md`
- `7_Treeshrink/dependencies/scripts/` 中保留了辅助脚本

### 8_1_concatenation

- 作用：拼接矩阵建树总流程
- 总说明：
- `8_1_concatenation/README.md`
- 步骤与脚本：
- `8_1_concatenation/0_concat/concat.py`：拼接多个 alignment
- `8_1_concatenation/1_del_indel/delmissingsite_all.py`：拼接矩阵缺失位点过滤
- `8_1_concatenation/2_codon/extract_and_merge_codons.py`：按密码子位置拆分和合并
- `8_1_concatenation/3_raxml/aln_no_outgroup.py`：RAxML 建树
- `8_1_concatenation/4_iqtree/run_iqtree_batch_simple.py`：IQ-TREE 建树
- `8_1_concatenation/5_mrbayes/run_mrbayes_workflow.py`：MrBayes 建树

### 8_2_coalescent

- 作用：共祖模型建树流程
- 总说明：
- `8_2_coalescent/README.md`
- 步骤与脚本：
- `8_2_coalescent/0_raxml/aln_no_outgroup.py`：逐基因建树
- `8_2_coalescent/1_filter/filter_all_in_one.py`：基因树过滤
- `8_2_coalescent/2_0_astral/run_astral.py`：ASTRAL 物种树推断

### 8_genetree_speciestree_conflict

- 作用：基因树与物种树冲突分析、定年配置生成、densitree 可视化
- 总说明：
- `8_genetree_speciestree_conflict/README.md`
- 步骤与脚本：
- `8_genetree_speciestree_conflict/0_extract_subtree/extract_subtrees.py`：提取目标子树
- `8_genetree_speciestree_conflict/1_generate_treePL_configfile/generate_treePL_conf.py`：生成 treePL 配置文件
- `8_genetree_speciestree_conflict/2_treePL/run_treepl.py`：运行 treePL
- `8_genetree_speciestree_conflict/3_draw_densitree/draw_densitree_workflow.py`：绘制 densitree 并统计节点支持

### 9_daing

- 作用：单次定年流程，围绕 treePL 展开
- 正式脚本：
- `9_daing/generate_treePL.py`
- 补充：
- `9_daing/README.md`

### 11_Hybridization_detection

- 作用：杂交检测总目录
- 总说明：
- `11_Hybridization_detection/README.md`
- 子目录：
- `11_Hybridization_detection/Phylonet/`：PhyloNet 路线
- `11_Hybridization_detection/QuIBL/`：QuIBL 路线

#### 11_Hybridization_detection/Phylonet

- 总说明：
- `11_Hybridization_detection/Phylonet/README.md`
- 步骤与脚本：
- `11_Hybridization_detection/Phylonet/0_将想要检测的树提取出来/extract_subtrees.py`
- `11_Hybridization_detection/Phylonet/1_Get_phylonet_input/newick2nex.py`
- `11_Hybridization_detection/Phylonet/2_Run_phylonet/aln.py`
- `11_Hybridization_detection/Phylonet/3_phylonet_result_to_dendroscope/get_tree_and_probability_from_phylonet.py`

#### 11_Hybridization_detection/QuIBL

- 总说明：
- `11_Hybridization_detection/QuIBL/README.md`
- 正式脚本：
- `11_Hybridization_detection/QuIBL/run_quibl.py`
- `11_Hybridization_detection/QuIBL/visual_quibl.py`

### 12_WGD_detection

- 作用：WGD 检测总目录
- 总说明：
- `12_WGD_detection/README.md`
- 子流程：
- `12_WGD_detection/absolute_dating/`：绝对定年路线
- `12_WGD_detection/phypart/`：PhyParts 路线

#### 12_WGD_detection/absolute_dating

- 总说明：
- `12_WGD_detection/absolute_dating/README.md`
- 步骤与脚本：
- `12_WGD_detection/absolute_dating/0_get_gene_family/get_gene_family.py`
- `12_WGD_detection/absolute_dating/1_mafft/mafft_pep_cds.py`
- `12_WGD_detection/absolute_dating/2_delmissing/delmissingsite.py`
- `12_WGD_detection/absolute_dating/3_raxml/aln_no_outgroup.py`
- `12_WGD_detection/absolute_dating/4_reroot_gene_family_trees/reroot_genefamily_by_species_tree.py`
- `12_WGD_detection/absolute_dating/5_extract_time/extract_time.py`
- `12_WGD_detection/absolute_dating/6_generate_treePL_file/generate_treePL_file.py`
- `12_WGD_detection/absolute_dating/7_visualize_results/run_result_plots.py`

#### 12_WGD_detection/phypart

- 总说明：
- `12_WGD_detection/phypart/README.md`
- 步骤与脚本：
- `12_WGD_detection/phypart/0_get_gene_family/get_gene_family.py`
- `12_WGD_detection/phypart/1_mafft/mafft_pep_cds.py`
- `12_WGD_detection/phypart/2_del_indel/delmissingsite.py`
- `12_WGD_detection/phypart/3_raxml/aln_no_outgroup.py`
- `12_WGD_detection/phypart/4_reroot_gene_family_trees/reroot_genefamily_by_species_tree.py`
- `12_WGD_detection/phypart/5_rename/rename.py`
- `12_WGD_detection/phypart/6_phypart/run_phypart.py`

### 20_Chloroplast

- 作用：叶绿体相关分析总目录
- 总说明：
- `20_Chloroplast/README.md`
- 子模块与入口脚本：
- `20_Chloroplast/0_密码子偏好/codon_usage_workflow.py`：密码子使用偏好分析
- `20_Chloroplast/1_ssr分析/ssr_workflow.py`：SSR 检测流程
- `20_Chloroplast/2_基因正选择分析/positive_selection_workflow.py`：叶绿体基因正选择分析
- `20_Chloroplast/3_基因表格生成/analyze_gb_genes.py`：基因组成和表格统计
- `20_Chloroplast/4_高变区/pi_workflow.py`：Pi 统计与高变区分析
- `20_Chloroplast/5_mvista/`：主要是 mVISTA 相关说明和示例文件，当前没有顶层 Python 入口脚本
- `20_Chloroplast/6_gb序列重命名/rename_locus.py`：GenBank 位点重命名
- `20_Chloroplast/7_chloroplast组装/chloroplast_assembly.py`：叶绿体组装流程
- `20_Chloroplast/8_叶绿体用cds建树/extract_orthologs.py`：提取同源 CDS 用于建树
- `20_Chloroplast/9_复杂基因结构建树/README.md`：复杂基因结构建树总流程说明
- `20_Chloroplast/PGA/PGA.pl`：保留的上游注释工具，不是本仓库重新封装后的 Python 主流程

#### 20_Chloroplast/9_复杂基因结构建树

- 步骤与脚本：
- `20_Chloroplast/9_复杂基因结构建树/0_把gb分割成小片段/extract_gene_region.py`
- `20_Chloroplast/9_复杂基因结构建树/1_proteinortho/run_proteinortho.py`
- `20_Chloroplast/9_复杂基因结构建树/2_split_tsv_by_type/split_tsv_by_type.py`
- `20_Chloroplast/9_复杂基因结构建树/3_parse_ortho_dna_groups/parse_ortho_dna_groups.py`
- `20_Chloroplast/9_复杂基因结构建树/4_mafft/mafft_and_concat.py`
- `20_Chloroplast/9_复杂基因结构建树/5_delmissing/delmissingsite.py`
- `20_Chloroplast/9_复杂基因结构建树/6_iqtree/run_iqtree.py`
- `20_Chloroplast/9_复杂基因结构建树/7_exabayes/run_exabayes.py`

### 21_annotation

- 作用：注释相关分析
- 总说明：
- `21_annotation/README.md`
- 子模块与入口脚本：
- `21_annotation/Genoma/run_gemoma.py`：GeMoMa 注释流程
- `21_annotation/annotate_fasta_genes/annotate_fasta_files.py`：给 FASTA 中的基因名称做注释转移或重命名

### 22_re-sequencing

- 作用：重测序相关分析总目录
- 总说明：
- `22_re-sequencing/README.md`
- 子流程：
- `22_re-sequencing/GetGeneFromGFF/`：从 GFF 与参考序列中提取目标基因
- `22_re-sequencing/ReSeqPipline/`：完整重测序分析主流程

#### 22_re-sequencing/GetGeneFromGFF

- 总说明：
- `22_re-sequencing/GetGeneFromGFF/README.md`
- 步骤与脚本：
- `22_re-sequencing/GetGeneFromGFF/0_filter_mRNA_annotation_from_gff/extract_mRNA_annotation.py`
- `22_re-sequencing/GetGeneFromGFF/1_get_gene_from_fasta_in_gff/gff2gene.py`
- `22_re-sequencing/GetGeneFromGFF/2_collection_gene/collection_genes.py`
- `22_re-sequencing/GetGeneFromGFF/3_del_missing/del_missing.py`

#### 22_re-sequencing/ReSeqPipline

- 总说明：
- `22_re-sequencing/ReSeqPipline/README.md`
- 主入口：
- `22_re-sequencing/ReSeqPipline/pipeline_master.py`
- 补充：
- `22_re-sequencing/ReSeqPipline/dependencies/scripts/` 中保留了分步骤脚本，如参考序列建索引、BWA 比对、GATK 调用、VCF 转 FASTA、序列合并等

### Other_tools

- 作用：通用辅助工具总目录
- 总说明：
- `Other_tools/README.md`
- 分类与脚本：

#### Other_tools/0_树相关

- `Other_tools/0_树相关/0_删掉树中的一些类群/del_tree.py`
- `Other_tools/0_树相关/0_删掉树中的一些类群/pure_tree.py`
- `Other_tools/0_树相关/1_提取树中的所有物种名称/aln.py`
- `Other_tools/0_树相关/2_重命名树中的所有样品名称/rename_newick_leaves.py`
- `Other_tools/0_树相关/3_对树重新置根/reroot.py`
- `Other_tools/0_树相关/4_输入物种名称随机生成进化树/generate_trees.py`
- `Other_tools/0_树相关/5_newick树上的bootstrap转为后验概率/process_tree.py`

#### Other_tools/1_统计各种数字

- `Other_tools/1_统计各种数字/0_fasta序列N50统计/assembly_stats.py`
- `Other_tools/1_统计各种数字/1_统计get_proteinortho结果中各个物种的基因数目/count_species.py`
- `Other_tools/1_统计各种数字/3_alignment整体信息位点统计/alignment_stats.py`
- `Other_tools/1_统计各种数字/4_alignment各序列gap信息统计/analyze_fasta.py`

#### Other_tools/2_画各种图

- `Other_tools/2_画各种图/0_频率分布直方图/plot_histogram.py`
- `Other_tools/2_画各种图/1_小提琴图/create_violin_plot.py`

#### Other_tools/3_fasta&fastq相关

- `Other_tools/3_fasta&fastq相关/0_重命名fasta中的所有样品名称/rename_fasta_headers.py`
- `Other_tools/3_fasta&fastq相关/1_bootstrap_fasta/bootstrap_generator.py`
- `Other_tools/3_fasta&fastq相关/2_提取fasta矩阵中的部分位点/extract_sites.py`
- `Other_tools/3_fasta&fastq相关/3_提取fasta矩阵中的部分序列/extract_by_name.py`
- `Other_tools/3_fasta&fastq相关/4_将一个大的fastq文件分割开/split_fastq.py`

#### Other_tools/4_arcgis相关

- `Other_tools/4_arcgis相关/0_掩膜提取/mask_and_warp.py`
- `Other_tools/4_arcgis相关/1_栅格图层格式转换/raster_converter.py`

#### Other_tools/5_sra数据提交

- `Other_tools/5_sra数据提交/ngdc/calculate_md5_from_list.py`

### lib

- 作用：仓库内部公用函数库
- 文件：
- `lib/__init__.py`
- `lib/utils.py`
- 说明：
- 这个目录主要供仓库中的脚本调用，不是面向最终用户的独立分析流程

## 快速定位指南

如果你已经知道自己要做什么，可以直接按下面方式进入：

- 从原始转录组数据开始跑主流程：`0_trimmomatic` 到 `7_Treeshrink`
- 做拼接树：`8_1_concatenation`
- 做共祖树：`8_2_coalescent`
- 做基因树/物种树冲突分析：`8_genetree_speciestree_conflict`
- 做单次定年：`9_daing`
- 做杂交检测：`11_Hybridization_detection`
- 做 WGD 检测：`12_WGD_detection`
- 做叶绿体分析：`20_Chloroplast`
- 做注释：`21_annotation`
- 做重测序：`22_re-sequencing`
- 找零散小工具：`Other_tools`

## 说明与边界

- 本文重点是“目录导航”和“脚本定位”，不展开每个脚本的参数细节
- 真正运行前，仍建议进入对应目录阅读本地 `README.md`
- `dependencies/` 中经常包含模板、上游说明、第三方程序或辅助脚本，通常不作为首选入口
- 某些目录保留了上游原始程序，例如 `3_TransDecoder-master/`、`20_Chloroplast/PGA/`，阅读时要区分“仓库封装脚本”和“上游工具本体”
