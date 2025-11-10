import argparse
import os


def index(reference):
    os.system("bwa index -a bwtsw " + reference)
    os.system("samtools faidx " + reference)

def main():
    #解析参数
    parser = argparse.ArgumentParser(description="Options for 1_Index_reference.py", add_help=True)
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-r', '--reference', action="store", metavar='\b', type=str, required=True, help="fasta格式的参考基因组")

    args                  = parser.parse_args()
    reference                = args.reference

    index(reference) #主程序   


main()