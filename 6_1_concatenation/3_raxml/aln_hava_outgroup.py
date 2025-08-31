from multiprocessing import Pool
import os, time, random

'''
定义参数
'''
n_cpu = 40   # 同时运行的进程数
Threads = 1  # 单个进程可以使用的线程数
file_tag = ".fas"  # 输入的文件(alignment)后缀
bootstraps = 1  # bootstraps次数
model = "GTRGAMMA"  # 使用的模型
outgroup = "Tupaia_chinensis_GCF_000334495.pep"  # 外类群的名称

def get_file_list():
    '''
    获取当前文件夹中符合目标扩展名的文件
    '''      
    file_name = []       
    for each in os.listdir(os.getcwd()):        
        if each[-len(file_tag):] == file_tag:
            file_name.append(each)
    return file_name


def get_th_list(file_name):
    '''
    get_th_list: 将所有文件平均的分给各个进程
    '''
    th_list = []
    for each_num in range(n_cpu):
        th_list.append([])
    print("本程序共使用 " + str(len(th_list)) + " 个进程")
    n = 0
    for each_file in file_name:       
        th_list[n].append(each_file)
        if n == n_cpu - 1:
            n = 0
        else:
            n = n + 1
    return th_list


def main_software(each_file_name_list):
    '''
    main_software: 运行raxml
    '''
    for each in each_file_name_list:
        command = ("raxmlHPC-PTHREADS" + 
                    " -T " + str(Threads) + 
                    " -n "+ each.replace(file_tag, "") +  
                    " -s " + each + 
                    " -m " + model + 
                    " -p 12345" + 
                    " -f a" + 
                    " -N " + str(bootstraps) +  
                    " -x 12345" +
                    " -o " + outgroup)  # 添加外类群参数
        print("run command " + command)
        os.system(command)


file_name = get_file_list()

th_list = get_th_list(file_name)

p = Pool(n_cpu)
for i in range(n_cpu):
    p.apply_async(main_software,(th_list[i],))

print("----start----")
p.close()  # 关闭进程池，关闭后po不再接收新的请求
p.join()  # 等待po中所有子进程执行完成
print("-----end-----")
