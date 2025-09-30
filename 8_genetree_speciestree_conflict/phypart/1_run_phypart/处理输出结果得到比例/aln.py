from ete3 import Tree


phypart_output_file = "out.concon.tre"
def main():
    with open(phypart_output_file, "r") as read_file:
        t1 = Tree(read_file.readline())
        t2 = Tree(read_file.readline())
        for each_node in t1.traverse():
            if each_node.is_leaf():
                pass
            else:
                name_list1 = each_node.get_leaf_names()
                name_list1.sort()
                for each_node2 in t2.traverse():
                    name_list2 = each_node2.get_leaf_names()
                    name_list2.sort()
                    if name_list1 == name_list2:
                        a = (each_node.support + each_node2.support)
                        if a == 0:
                            a = 0.00001
                        each_node.support = each_node.support/a
                        #print(str(each_node.support) + " " + str(each_node.support/(each_node.support + each_node2.support)))
    t1.write(format=0, outfile="new_tree.nw")





main()