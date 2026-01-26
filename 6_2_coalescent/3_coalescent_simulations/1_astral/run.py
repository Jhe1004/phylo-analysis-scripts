import os
from ete3 import Tree

#os.system("java -jar astral.jar -i output_subtree.trees -o astral.tree")

tree1 = Tree("astral.tree")

tree1.write(format=1, outfile="astral.tre")

os.system("perl unifyTreeLength.pl astral.tre")

os.remove("astral.tree")

os.remove("astral.tre")