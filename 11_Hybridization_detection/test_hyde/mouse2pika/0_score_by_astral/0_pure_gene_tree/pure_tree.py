import os
from ete3 import Tree
now_dir = os.getcwd()    

ingroups_file = "pika.tree.txt"

def get_file_list():
    file_temp = os.listdir()
    file_list = []
    for each in file_temp:       
        if "RAx" in each:
            file_list.append(each)
    return file_list


def main(inpurt_tree):
    ingroups_list = []
    with open(ingroups_file, "r") as read_file:
        for each_line in read_file:
            ingroups_list.append(each_line.replace("\n", ""))
    t = Tree(inpurt_tree)
    name_list_in_genetree = t.get_leaf_names()
    name_list_in_genetree_new = []
    for each_name in name_list_in_genetree:
        if each_name in ingroups_list:
            name_list_in_genetree_new.append(each_name)
        else:
            pass
    t.prune(name_list_in_genetree_new)
    t.write(format=1, outfile=inpurt_tree + "_new.tree")
    for each_name in ingroups_list:
        if each_name not in name_list_in_genetree_new:
            print(each_name + "not in the subtree")


if __name__ == "__main__":
    for each_tree_file in get_file_list():
        main(each_tree_file)