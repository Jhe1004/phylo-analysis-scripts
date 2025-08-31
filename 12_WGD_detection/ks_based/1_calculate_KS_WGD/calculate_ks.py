'''
本脚本并行计算所有物种的ks值。
本脚本使用python3运行，但脚本中的功能需要调用python2环境，并且该python2环境中需要安装biopython
'''

from multiprocessing import Pool
import os
import subprocess

'''
首先定义需要的进程数th
'''
th = 35


'''
函数 get_file_list: 获取当前文件夹中符合目标扩展名的文件
输入 无，将本脚本放置在目标文件夹中即可
输出 file_name：所有文件的名称列表
'''
def get_file_list():      
    file_name = []       
    for each in os.listdir(os.getcwd()):        
        if ".cds" in each:
            file_name.append(each.replace(".cds", ""))
    return file_name

'''
函数 get_th_list: 将所有文件平均的分给各个进程
输入 file_name: 所有文件的名称列表
输出 th_list： 一个列表，列表中的各个元素为各个进程所应该处理的文件
'''
def get_th_list(file_name):
    th_list = []
    for each_num in range(th):
        th_list.append([])
    print("本程序共使用 " + str(len(th_list)) + " 个进程")
    n = 0
    for each_file in file_name:       
        th_list[n].append(each_file)
        if n == th - 1:
            n = 0
        else:
            n = n + 1
    return th_list


def main_software(each_file_name_list):
    print(each_file_name_list)
    for each in each_file_name_list:
        subprocess.call("mkdir " + each, shell=True)
        subprocess.call("cp " + each + "* ./" + each, shell=True)
        subprocess.call("cp kSPlotter.py " + "./" + each, shell=True)
        work_dir =  os.getcwd()
        os.chdir(work_dir + "/" + each)
        subprocess.call("python kSPlotter.py -b " + each + ".pepresult -aa " + each + ".pep -nt " + each + ".cds -o " + each, shell=True)
        os.chdir(work_dir)        

'''
#无法运行时单进程调试代码,注释掉本框后所有部分使用
file_name = get_file_list()
th_list = get_th_list(file_name)
for each in th_list:
    main_software(each)

'''
file_name = get_file_list()
th_list = get_th_list(file_name)

p = Pool(th)
for i in range(th):
    p.apply_async(main_software,(th_list[i],))

print("----start----")
p.close()  # 关闭进程池，关闭后po不再接收新的请求
p.join()  # 等待po中所有子进程执行完成，再执行下面的代码,可以设置超时时间join(timeout=)
print("-----end-----") 
