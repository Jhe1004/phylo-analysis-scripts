import os
import argparse  


'''
获取文件夹中的样品名称
'''


def main():
    #解析参数
    parser = argparse.ArgumentParser(description="Options for 4_get_fasta.py", add_help=True)
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-r', '--reference', action="store", metavar='\b', type=str, required=True, help="fasta格式的参考基因组")
    additional = parser.add_argument_group("Additional arguments")
    additional.add_argument('-f', '--folder', action="store", metavar='\b', type=str,  default="current_location", help="存放vcf文件的文件夹路径，默认为本脚本所在的文件夹")
    
    
    args                = parser.parse_args()
    reference           = args.reference
    folder              = args.folder



    if folder == "current_location":
        folder = os.getcwd()
    file_tag = ".vcf"
    res_list = []
    for each_file in os.listdir(folder):       
        if file_tag == each_file[-len(file_tag):]:
            res_list.append(each_file)
    
    with open(reference, "r") as read_file:
        for each_line in read_file:
            if each_line[0] == ">":
                contig_name = each_line.split(" ")[0][1:].replace("\n","")
                with open(contig_name + "_combine.fasta", "a") as write_file:
                    for each_sample in res_list:
                        with open(contig_name + "_" + each_sample + ".fasta", "r") as read_file2:
                            for each_line in read_file2:
                                write_file.write(each_line)
   
    
    

main()
