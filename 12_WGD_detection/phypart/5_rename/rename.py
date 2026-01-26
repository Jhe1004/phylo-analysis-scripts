import os
import re

'''
函数 get_file_list: 获取当前文件夹中符合目标扩展名的文件
输入 无，将本脚本放置在目标文件夹中即可
输出 file_name：所有文件的名称列表
'''

file_tag = ".tree"
def get_file_list():
    '''
    获得当前文件夹中相应后缀的文件列表
    '''   
    res_list = []
    for each_file in os.listdir():       
        if file_tag == each_file[-len(file_tag):]:
            res_list.append(each_file)
    return res_list

file_name = get_file_list()
for each in file_name:
    with open(each.replace(file_tag, "_re.tree"), "a") as write_file:
        with open(each, "r") as read_file:
            for each_line in read_file:
                    pattern1 = "\+\+.+?\:"
                    newseq = re.sub(pattern1, ":", each_line)
                    write_file.write(newseq)                    


