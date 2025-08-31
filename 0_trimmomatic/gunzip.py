from multiprocessing import Pool
import os, time, random

'''
首先定义需要的进程数th
'''
th = 39


'''
函数 get_file_list: 获取当前文件夹中符合目标扩展名的文件
输入 无，将本脚本放置在目标文件夹中即可
输出 file_name：所有文件的名称列表
'''
def get_file_list():      
    file_name = []       
    for each in os.listdir(os.getcwd()):        
        if ".gz" in each:
            file_name.append(each)
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
        os.system("gunzip "+ each)



file_name = get_file_list()
th_list = get_th_list(file_name)

p = Pool(th)
for i in range(th):
    p.apply_async(main_software,(th_list[i],))

print("----start----")
p.close()  # 关闭进程池，关闭后po不再接收新的请求
p.join()  # 等待po中所有子进程执行完成，再执行下面的代码,可以设置超时时间join(timeout=)
print("-----end-----") 
  
		
           
