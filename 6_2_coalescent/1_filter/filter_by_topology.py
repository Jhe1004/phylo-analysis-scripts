from ete3 import Tree
import os

# 定义物种树的文件路径或直接使用Newick字符串
species_tree_newick = "((Aletris_farinosa_SRR19179843.fasta.transdecoder.pep:1.0,(Burmannia_biflora_SRR19179841.fasta.transdecoder.pep:1.0,Dioscorea_oppositifolia_SRR25002067.fasta.transdecoder.pep:1.0)100.0:1.0)100.0:0.5,((Acanthochlamys_bracteata_GCA_019914995.fasta.transdecoder.pep:1.0,(Xerophyta_villosa_ERR2040709.fasta.transdecoder.pep:1.0,(Xerophyta_humilis_SRR8334914.fasta.transdecoder.pep:1.0,Talbotia_elegans_SRR19179852.fasta.transdecoder.pep:1.0)100.0:1.0)100.0:1.0)100.0:1.0,((Lacandonia_schismatica_SRR19179849.fasta.transdecoder.pep:1.0,(Triuris_sp_2_SRR19179851.fasta.transdecoder.pep:1.0,Triuris_sp_1_SRR19179850.fasta.transdecoder.pep:1.0)100.0:1.0)100.0:1.0,((Stemona_collinsiae_SRR15465044.fasta.transdecoder.pep:1.0,(Croomia_pauciflora_SRR8298337.fasta.transdecoder.pep:1.0,(Croomia_heterosepala_SRR6220063.fasta.transdecoder.pep:1.0,Croomia_japonica_SRR6219921.fasta.transdecoder.pep:1.0)100.0:1.0)100.0:1.0)100.0:1.0,(Sphaeradenia_woodsonii_ERR3487374.fasta.transdecoder.pep:1.0,(Freycinetia_multiflora_ERR2040706.fasta.transdecoder.pep:1.0,(Pandanus_amaryllifolius_SRR13452149.fasta.transdecoder.pep:1.0,(Pandanus_utilis_SRR6374687.fasta.transdecoder.pep:1.0,(Pandanus_odorifer_SRR5927127.fasta.transdecoder.pep:1.0,Pandanus_tectorius_SRR19179847.fasta.transdecoder.pep:1.0)100.0:1.0)100.0:1.0)100.0:1.0)100.0:1.0)100.0:1.0)100.0:1.0)100.0:1.0):0.5);"

# 读取物种树
species_tree = Tree(species_tree_newick, format=1)

# 定义基因树所在的目录
gene_trees_dir = "./"  # 请替换为实际路径

# 定义保留差异阈值（Robinson-Foulds距离）
rf_distance_threshold = 10  # 可以根据需要调整

# 创建一个用于保存保留基因树的目录
output_dir = "filtered_gene_trees"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 遍历基因树文件
for filename in os.listdir(gene_trees_dir):
    if filename.endswith(".nwk") or filename.endswith(".newick") or filename.endswith(".tree"):
        gene_tree_path = os.path.join(gene_trees_dir, filename)
        gene_tree = Tree(gene_tree_path, format=1)
        
        # 获取物种树和基因树的叶节点名称集合
        species_leaves = set(species_tree.get_leaf_names())
        gene_leaves = set(gene_tree.get_leaf_names())

        # 获取公共的叶节点名称
        common_leaves = species_leaves & gene_leaves

        # 如果公共叶节点数量不足，可以选择跳过该基因树
        if len(common_leaves) < 3:  # 最少需要3个公共叶节点
            print(f"{filename} 跳过，公共叶节点数量不足")
            continue

        # 创建物种树和基因树的副本，以免修改原始树
        species_tree_copy = species_tree.copy()
        gene_tree_copy = gene_tree.copy()

        # 修剪树，只保留公共叶节点
        species_tree_copy.prune(common_leaves, preserve_branch_length=True)
        gene_tree_copy.prune(common_leaves, preserve_branch_length=True)

        # 计算RF距离
        rf_result = species_tree_copy.robinson_foulds(gene_tree_copy, unrooted_trees=True)
        rf_distance = rf_result[0]
        max_rf = rf_result[1]

        # 如果距离小于阈值，则保留基因树
        if rf_distance <= rf_distance_threshold:
            output_path = os.path.join(output_dir, filename)
            gene_tree.write(format=1, outfile=output_path)
            print(f"{filename} 保留，RF距离：{rf_distance}")
        else:
            print(f"{filename} 丢弃，RF距离：{rf_distance}")