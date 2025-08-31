import os


simulate_times  = 1000  #模拟多少次
species_num = 20 #物种树中的样品数目
species_tree = "coalescent_units_tree.newick"

    
#使用R脚本进行模拟
os.system("Rscript 00.simulate.Rscript " + species_tree + " " + str(species_num) + " " + str(simulate_times))


