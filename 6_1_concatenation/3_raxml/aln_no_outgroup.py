from multiprocessing import Pool
import os

'''
定义参数
'''
n_cpu = 1           # 同时运行的进程数
Threads = 10        # 单个进程可以使用的线程数
file_tag = ".fasta" # 输入的文件后缀
bootstraps = 10     # bootstraps次数
model = "GTRGAMMA"  # 使用的模型

def get_file_list():
    '''获取当前文件夹中符合目标扩展名的文件'''      
    file_name = []       
    for each in os.listdir(os.getcwd()):        
        if each.endswith(file_tag): # 使用 endswith 更稳健
            file_name.append(each)
    return file_name

def get_th_list(file_name, actual_n_cpu):
    '''将所有文件平均的分给各个进程'''
    th_list = [[] for _ in range(actual_n_cpu)]
    
    n = 0
    for each_file in file_name:       
        th_list[n].append(each_file)
        n = (n + 1) % actual_n_cpu # 自动循环索引
    return th_list

def main_software(each_file_name_list):
    '''运行 raxml'''
    for each in each_file_name_list:
        # 移除后缀获取前缀
        prefix = each.replace(file_tag, "")
        command = (f"raxmlHPC -T {Threads} -n {prefix} "
                   f"-s {each} -m {model} -p 12345 -f a "
                   f"-N {bootstraps} -x 12345")
        print(f"Running command: {command}")
        os.system(command)

# --- 关键修改：必须放在这个判断下面 ---
if __name__ == '__main__':
    # 1. 获取文件列表
    file_list = get_file_list()
    
    if not file_list:
        print(f"未找到后缀为 {file_tag} 的文件，请检查目录。")
    else:
        # 2. 确定实际使用的进程数（不能超过文件总数）
        actual_n_cpu = min(len(file_list), n_cpu)
        print(f"本程序共使用 {actual_n_cpu} 个进程")

        # 3. 分配任务
        th_list = get_th_list(file_list, actual_n_cpu)

        # 4. 启动进程池
        p = Pool(actual_n_cpu)
        for i in range(actual_n_cpu):
            p.apply_async(main_software, (th_list[i],))

        print("---- Start multiprocessing ----")
        p.close()  # 关闭进程池，不再接收新任务
        p.join()   # 等待所有子进程完成
        print("---- All tasks finished ----")