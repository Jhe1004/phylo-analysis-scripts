import os


simulate_times  = 10000  #模拟多少次
species_num = 27 #物种树中的样品数目
species_tree = "astral.tre.regular.tre"

    
#使用R脚本进行模拟
os.system("Rscript 00.simulate.Rscript " + species_tree + " " + str(species_num) + " " + str(simulate_times))


