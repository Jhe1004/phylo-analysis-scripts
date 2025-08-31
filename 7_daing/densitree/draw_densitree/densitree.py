import toytree
import toyplot
import numpy as np
from matplotlib import pyplot as plt
import toyplot.svg


# set dimensions of the canvas
canvas = toyplot.Canvas(width=2400, height=750)


ax0 = canvas.cartesian(bounds=('15%', '44.6%', '10%', '90%'))
ax1 = canvas.cartesian(bounds=('15%', '44.6%', '10%', '90%'))


mtrees = toytree.mtree("result.trees") #基因树
tree1 = toytree.tree("dating.tree") #物种树
species_order = ['Anemoclema_glauciifolium_19.fasta.transdecoder.pep', 'Clematis_huchouensis.fasta.transdecoder.pep', 'Clematis_fusca.fasta.transdecoder.pep', 'Clematis_reticulata.fasta.transdecoder.pep', 'Clematis_loureiriana.fasta.transdecoder.pep', 'Clematis_uncinata.fasta.transdecoder.pep', 'Clematis_lancifolia.fasta.transdecoder.pep', 'Clematis_armandii.fasta.transdecoder.pep', 'Clematis_crassifolia.fasta.transdecoder.pep', 'Clematis_quinquefoliolata.fasta.transdecoder.pep', 'Clematis_terniflora.fasta.transdecoder.pep', 'Clematis_hexapetala.fasta.transdecoder.pep', 'Clematis_macropetala.fasta.transdecoder.pep', 'Clematis_ochotensis_2019051801.fasta.transdecoder.pep','Clematis_acerifolia.fasta.transdecoder.pep', 'Clematis_aethusifolia.fasta.transdecoder.pep', 'Clematis_peterae.fasta.transdecoder.pep', 'Clematis_gouriana.fasta.transdecoder.pep', 'Clematis_gratopsis.fasta.transdecoder.pep', 'Clematis_brevicaudata.fasta.transdecoder.pep', 'Clematis_tubulosa.fasta.transdecoder.pep', 'Clematis_puberula.fasta.transdecoder.pep', 'Clematis_delevayi.fasta.transdecoder.pep', 'Clematis_lasiandra.fasta.transdecoder.pep', 'Clematis_leschenaultiana2.fasta.transdecoder.pep', 'Clematis_leschenaultiana.fasta.transdecoder.pep', 'Clematis_rehderiana.fasta.transdecoder.pep', 'Naravelia_siamensis.fasta.transdecoder.pep', 'Clematis_montana.fasta.transdecoder.pep', 'Clematis_repens.fasta.transdecoder.pep', 'Clematis_otophora.fasta.transdecoder.pep', 'Clematis_pseudopogonandra.fasta.transdecoder.pep', 'Clematis_orientalis_kashi.fasta.transdecoder.pep', 'Clematis_glauca_liancheng.fasta.transdecoder.pep', 'Clematis_intricata_2019052001.fasta.transdecoder.pep', 'Clematis_intricata_weijing.fasta.transdecoder.pep', 'Clematis_tenuifolia.fasta.transdecoder.pep', 'Clematis_akebioides_ganzi.fasta.transdecoder.pep', 'Clematis_tangutica_yongchang.fasta.transdecoder.pep', 'Clematis_tangutica_maerkang.fasta.transdecoder.pep', 'Clematis_akebioides_wenchuan.fasta.transdecoder.pep', 'Clematis_viridis_mupo.fasta.transdecoder.pep', 'Clematis_viridis_jiangda.fasta.transdecoder.pep', 'Clematis_viridis_yushu.fasta.transdecoder.pep', 'Clematis_nannophylla_tanchang.fasta.transdecoder.pep', 'Clematis_nannophylla_zhouqu1.fasta.transdecoder.pep', 'Clematis_nannophylla_zhouqu2.fasta.transdecoder.pep', 'Clematis_nannophylla_wudu.fasta.transdecoder.pep', 'Clematis_nannophylla_wudu2.fasta.transdecoder.pep', 'Clematis_nannophylla_jianzha.fasta.transdecoder.pep', 'Clematis_nannophylla_mengda.fasta.transdecoder.pep', 'Clematis_nannophylla_baiyin.fasta.transdecoder.pep', 'Clematis_nannophylla_lanzhou.fasta.transdecoder.pep', 'Clematis_tomentella_zhangye.fasta.transdecoder.pep', 'Clematis_tomentella_jiayuguan.fasta.transdecoder.pep', 'Clematis_tomentella_lingwu.fasta.transdecoder.pep', 'Clematis_canescens_weijing.fasta.transdecoder.pep', 'Clematis_canescens_weijing2.fasta.transdecoder.pep', 'Clematis_fruticosa_baotou.fasta.transdecoder.pep', 'Clematis_fruticosa_yanan.fasta.transdecoder.pep', 'Clematis_fruticosa_luliang.fasta.transdecoder.pep', 'Clematis_fruticosa_datong.fasta.transdecoder.pep', 'Clematis_fruticosa_wutaishan.fasta.transdecoder.pep', 'Clematis_songorica_hami.fasta.transdecoder.pep', 'Clematis_songorica_kashi.fasta.transdecoder.pep', 'Clematis_songorica_wulumuqi.fasta.transdecoder.pep', 'Clematis_songorica_aleotai.fasta.transdecoder.pep']
species_order.reverse()
print(species_order)

mtrees.draw_cloud_tree(axes=ax0,
    edge_style={
        "stroke-opacity": 0.09,
        "stroke-width": 1.0,
    }, edge_type='c', edge_colors=toytree.colors[2], width=400, height=400, fixed_order=species_order, jitter=0.0)

tree1.draw(axes=ax1, edge_type='c', layout='r', edge_widths=0.3, fixed_order=species_order)

toyplot.svg.render(canvas, "tree-plot.svg")


