#!/usr/bin/python3
'''
运行bwa软件，将WGS序列map到基因组上

软件安装：
conda install bwa

'''

import os
import argparse




def map(reference, threads, plus_tag, minus_tag, folder):
    # 脚本需要批量运行，获得原始序列文件夹中所有的正向reads
    if folder == "current_location":
        folder = os.getcwd()
    file_name = []
    for each in os.listdir(folder):    
        if plus_tag == each[-(len(plus_tag)):]:
            file_name.append(each)
    # 运行bwa
    for each_fq in file_name:
        sample_name = each_fq.replace(plus_tag, "")
        plus_fq = folder + os.sep + each_fq
        minus_fq = folder + os.sep + each_fq.replace(plus_tag, minus_tag)
        map_command = ('''bwa mem -t %s -M -R \'@RG\\tID:%s\\tSM:%s\\tLB:WES\\tPL:Illumina\' %s %s %s > %s/%s.pre_sort.bam''' %(threads, sample_name, sample_name, reference, plus_fq, minus_fq, folder,sample_name))
        print("Map to Reference with bwa command: " + map_command)
        os.system(map_command)
        #对bam文件进行排序
        samtools_command2 = ('''samtools sort %s/%s.pre_sort.bam -o %s/%s.bam''' %(folder, sample_name, folder, sample_name))
        print("Sort BAM File with samtools command: " + samtools_command2)
        os.system(samtools_command2)
        #建立索引
        samtools_command3 = ('''samtools index -c %s/%s.bam''' %(folder, sample_name))
        print("Index BAM File with samtools command: " + samtools_command3)
        os.system(samtools_command3)
        os.system("rm *.pre_sort.bam")



def main():
    #解析参数
    parser = argparse.ArgumentParser(description="Options for 2_bwa_map.py", add_help=True)
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-r', '--reference', action="store", metavar='\b', type=str, required=True, help="fasta格式的参考基因组")
    required.add_argument('-t', '--threads', action="store", metavar='\b', type=str, required=True, help="程序最多可以使用的线程数目")  
    required.add_argument('-p', '--plus_tag', action="store", metavar='\b', type=str, required=True, help="正向测序数据的后缀，例如双端测序返回的原始数据名称分别为ginkgo_biloba_1.fastq.gz和ginkgo_biloba_2.fastq.gz，这里输入_1.fastq.gz即可")    
    required.add_argument('-m', '--minus_tag', action="store", metavar='\b', type=str, required=True, help="负向测序数据的后缀，例如双端测序返回的原始数据名称分别为ginkgo_biloba_1.fastq.gz和ginkgo_biloba_2.fastq.gz，这里输入_2.fastq.gz即可")     
    additional = parser.add_argument_group("Additional arguments")
    additional.add_argument('-f', '--folder', action="store", metavar='\b', type=str,  default="current_location", help="存放数据的文件夹路径，默认为本脚本所在的文件夹")
    
    
    args                = parser.parse_args()
    reference           = args.reference
    threads             = args.threads
    plus_tag            = args.plus_tag
    minus_tag           = args.minus_tag
    folder              = args.folder


    map(reference, threads, plus_tag, minus_tag, folder) #主程序   



main()