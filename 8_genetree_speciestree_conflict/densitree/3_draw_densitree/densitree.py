

import toytree
'''
绘制云图，云图中的实线为dating得到的物种树分化时间图
          云图中的虚线为dating得到的基因树分化时间图
需要输入物种树文件名以及基因树文件名
'''
import toyplot
import numpy as np
from matplotlib import pyplot as plt
import toyplot.svg


# set dimensions of the canvas
canvas = toyplot.Canvas(width=1900, height=1250)


ax0 = canvas.cartesian(bounds=('15%', '44.6%', '10%', '90%'))
ax1 = canvas.cartesian(bounds=('15%', '44.6%', '10%', '90%'))


mtrees = toytree.mtree("gene.trees")
tree1 = toytree.tree("dna.tree")
species_order = ["6_Calathodes_oxycarpa_2024050902-8_trinity.Trinity.fasta.transdecoder.pep", "A.r.fupingensis2_2023051702_trinity.Trinity.fasta.transdecoder.pep", "A.r.fupingensis4_2023051704_trinity.Trinity.fasta.transdecoder.pep", "A.r.fupingensis5_2023051705_trinity.Trinity.fasta.transdecoder.pep", "A_s_I-4749_I-4749_trinity.Trinity.fasta.transdecoder.pep", "A_sibirica_Altay_20230609035-6_trinity.Trinity.fasta.transdecoder.pep", "A_sibirica_4786_I-4786_trinity.Trinity.fasta.transdecoder.pep", "A_c_1_2024060508-1_trinity.Trinity.fasta.transdecoder.pep", "A_c_2_2024060508-2_trinity.Trinity.fasta.transdecoder.pep", "A_coerulea_2024060508_trinity.Trinity.fasta.transdecoder.pep", "A_sutchuenensis_1_WH25032502_trinity.Trinity.fasta.transdecoder.pep", "A_sutchuenensis_2_WH25032502-2_trinity.Trinity.fasta.transdecoder.pep", "Adonis_davidii_2025050306-12_trinity.Trinity.fasta.transdecoder.pep", "Adonis_davidii_2025050306-9_trinity.Trinity.fasta.transdecoder.pep", "Adonis_davidii_2025050306-1_trinity.Trinity.fasta.transdecoder.pep", "A_a_Tonghua_1_2024050405-1_trinity.Trinity.fasta.transdecoder.pep", "A_a_Tonghua_3_2024050405-3_trinity.Trinity.fasta.transdecoder.pep", "A_amurensis_Tonghua_2024050305-1_trinity.Trinity.fasta.transdecoder.pep", "A_amurensis_ChengnanPark2025031401_trinity.Trinity.fasta.transdecoder.pep", "A_amurensis_Japan_WHJP01_trinity.Trinity.fasta.transdecoder.pep", "A_r_2024050304-2_trinity.Trinity.fasta.transdecoder.pep", "A_r_Tonghua_1_2024050304-1_trinity.Trinity.fasta.transdecoder.pep", "A_r_Tonghua_5_2024050304-5_trinity.Trinity.fasta.transdecoder.pep", "A_r_SHandong_1_2025041001-1_trinity.Trinity.fasta.transdecoder.pep", "A_r_SHandong_5_2025041001-5_trinity.Trinity.fasta.transdecoder.pep", "A_r_Shandong-9_2025041001-9_trinity.Trinity.fasta.transdecoder.pep", "A_r_Henan-4_2025031201-4_trinity.Trinity.fasta.transdecoder.pep", "A_r_Henan_5_2025031201-5_trinity.Trinity.fasta.transdecoder.pep", "A_r_Henan_17_2025031201-17_trinity.Trinity.fasta.transdecoder.pep", "A_r_Huairou_Fenghuangtuo_2025032802-3_trinity.Trinity.fasta.transdecoder.pep", "A_r_chashikou2024030801_trinity.Trinity.fasta.transdecoder.pep", "A_r_Jiangsu_01_2025032301-1_trinity.Trinity.fasta.transdecoder.pep", "A_r_Jiangsu_03_2025032303-1_trinity.Trinity.fasta.transdecoder.pep", "A_bobro_1_2024060601-1_trinity.Trinity.fasta.transdecoder.pep", "A_bobro_2024060601_trinity.Trinity.fasta.transdecoder.pep", "A_bobro_3_2024060601-3_trinity.Trinity.fasta.transdecoder.pep", "A_sibirica_20180525_01.fasta.transdecoder.pep", "A_villosa_lll20240520011_trinity.Trinity.fasta.transdecoder.pep", "A_villosa_4247_I-4247_trinity.Trinity.fasta.transdecoder.pep", "A_aestivalis_var_parviflora_I-4422_trinity.Trinity.fasta.transdecoder.pep"]
#species_order = tree1.get_tip_labels()
print(species_order)

mtrees.draw_cloud_tree(axes=ax0,
    edge_style={
        "stroke-opacity": 0.09,
        "stroke-width": 1.0,
    }, edge_colors=toytree.colors[2], width=400, height=600, fixed_order=species_order)

tree1.draw(axes=ax1, edge_type='c', layout='r', edge_widths=2, fixed_order=species_order)

toyplot.svg.render(canvas, "tree-plot.svg")


