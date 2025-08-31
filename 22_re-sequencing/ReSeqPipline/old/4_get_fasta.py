from multiprocessing import Pool
import os
import argparse
import subprocess


#Divide all files to processes 
def get_th_list(file_name, threads):
    th_list = []
    for each_num in range(threads):
        th_list.append([])
    print("本程序共使用 " + str(len(th_list)) + " 个进程")
    n = 0
    for each_file in file_name:       
        th_list[n].append(each_file)
        if n == threads - 1:
            n = 0
        else:
            n = n + 1
    return th_list


def main_software(each_file_name_list,reference,folder):
    print(each_file_name_list)
    for each_vcf in each_file_name_list:
        command = ("python gvcf2fasta.py -r " + reference + " -g " + (folder + os.sep + each_vcf))
        print("Run gvcf2fasta with command " + command)
        subprocess.call(command, shell=True)

        




def main():
    #解析参数
    parser = argparse.ArgumentParser(description="Options for 3_gatk.py", add_help=True)
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-r', '--reference', action="store", metavar='\b', type=str, required=True, help="fasta格式的参考基因组")
    required.add_argument('-t', '--threads', action="store", metavar='\b', type=int, required=True, help="程序最多可以使用的线程数目")  
    additional = parser.add_argument_group("Additional arguments")
    additional.add_argument('-f', '--folder', action="store", metavar='\b', type=str,  default="current_location", help="存放vcf文件的文件夹路径，默认为本脚本所在的文件夹")
    
    
    args                = parser.parse_args()
    reference           = args.reference
    threads             = args.threads
    folder              = args.folder

    if folder == "current_location":
        folder = os.getcwd()
    file_name = []       
    for each in os.listdir(folder):        
        if each[-3:] == "vcf":
            file_name.append(each)

    th_list = get_th_list(file_name, threads)

    p = Pool(threads)
    for i in range(threads):
        p.apply_async(main_software,(th_list[i],reference,folder))

    print("----start----")
    p.close()  
    p.join()  
    print("-----end-----") 
    
    


main()

    