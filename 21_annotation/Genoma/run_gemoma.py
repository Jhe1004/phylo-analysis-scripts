'''
运行GeMoMa软件对基因组进行注释

软件官网：http://www.jstacs.de/index.php/GeMoMa

软件安装：
conda create -n gemoma
conda activate gemoma
conda install -c bioconda gemoma

软件输入：
需要注释的基因组：target_genome
参考基因组的fasta序列：reference_1_genome
参考基因组的gff格式注释文件：reference_1_annotation
运行命令为：
GeMoMa -Xmx100G GeMoMaPipeline threads=<threads> outdir=<outdir> GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=<target_genome> i=<reference_1_id> a=<reference_1_annotation> g=<reference_1_genome>
注：java默认的内存一般不够用，所以要使用参数-Xmx100G来指定java可以使用的最高内存大小
'''

import os
import glob

theads = "70"
Xmx = "50G"
reference_1_genome = "GCF_000001635.27_GRCm39_genomic.fasta"
reference_1_annotation = "genomic.gff"

for each_target_genome in glob.glob("*.fna"):
    outdir = each_target_genome.split(".fna")[0]
    target_genome = each_target_genome
    print(outdir)
    print(target_genome)

    command = ("GeMoMa -Xmx%s GeMoMaPipeline threads=%s outdir=%s GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true pc=true t=%s a=%s g=%s" %(Xmx,theads,outdir,target_genome,reference_1_annotation,reference_1_genome))


          
    print("run command" + command)
    os.system(command)
 
