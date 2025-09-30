import os
from Bio import SeqIO
from Bio.SeqUtils import GC
import re


file_tag = ".gb"  # Specify the file extension to process

# Function to get a list of files with the specified extension
def get_file_list():
    '''
    Get a list of files in the current directory with the specified extension.
    '''   
    res_list = []
    for each_file in os.listdir():       
        if each_file.endswith(file_tag):
            res_list.append(each_file)
    return res_list

def find_ir_region(seq):
    minimum_repeat_length = 1000
    irb_length = minimum_repeat_length
    judge = 0
    while judge != -1:
        irb_length += 1
        ira_left_site = judge
        # Find the reverse complement of the current IR region
        search_seq = seq[-irb_length:].reverse_complement()
        judge = seq.find(search_seq)
    irb_right_site = len(seq)
    irb_left_site = irb_right_site - irb_length + 2
    ira_right_site = ira_left_site + irb_length - 1
    ira_left_site = ira_left_site + 1
    return irb_right_site, irb_left_site, ira_right_site, ira_left_site             

def get_assembly_name(filename):
    with open(filename, 'r') as file:
        # 读取文件内容
        content = file.read()

    # 查找 Assembly Name 在 Assembly Data块中的位置
    match = re.search(r"##Assembly-Data-START##(.*?)##Assembly-Data-END##", content, re.DOTALL)
    
    if match:
        # 提取出 Assembly-Data 块的内容
        assembly_data_block = match.group(1)
        
        # 查找 Assembly Name 行
        assembly_name_match = re.search(r"Assembly Name\s*::\s*(\S+)", assembly_data_block)
        
        if assembly_name_match:
            return assembly_name_match.group(1)
        else:
            return "unknown"
    else:
        return "unknown"

def get_coding_size(record):
    temp_seq = "a" * len(record.seq)  # Temporary sequence to mark CDS regions
    
    for feature in record.features:
        if feature.type == "CDS":  # Only handle CDS features
            for part in feature.location.parts:  # Check each part of the location
                # Mark the corresponding regions as 'b' (coding region)
                left, right = part.start, part.end
                temp_seq = temp_seq[:left] + "b" * (right - left) + temp_seq[right:]

    # Count the number of 'b' characters, which represent the CDS regions
    return temp_seq.count("b")


def get_gene_num(record, ira_right_site, ira_left_site):
    gene_list = []
    cds_list = []
    trna_list = []
    rrna_list = []
    ir_gene_list = []
    
    for each_annotation in record.features:
        # Handle 'gene' features
        if each_annotation.type == "gene":
            gene_names = each_annotation.qualifiers.get("gene", [])
            for gene in gene_names:
                if gene not in gene_list:
                    gene_list.append(gene)
        
        # Handle 'CDS' features
        if each_annotation.type == "CDS":
            gene_names = each_annotation.qualifiers.get("gene", [])
            for gene in gene_names:
                if gene not in cds_list:
                    cds_list.append(gene)
        
        # Handle 'tRNA' features
        if each_annotation.type == "tRNA":
            gene_names = each_annotation.qualifiers.get("gene", [])
            for gene in gene_names:
                if gene not in trna_list:
                    trna_list.append(gene)
        
        # Handle 'rRNA' features
        if each_annotation.type == "rRNA":
            gene_names = each_annotation.qualifiers.get("gene", [])
            for gene in gene_names:
                if gene not in rrna_list:
                    rrna_list.append(gene)
        
        # Handle genes located within the inverted repeat (IR) regions
        if each_annotation.type == "gene":
            temp = str(each_annotation.location).replace("[","").replace("](+)","").replace("](-)","")
            gene_names = each_annotation.qualifiers.get("gene", [])
            if "join" not in temp:
                try:
                    left_site = int(temp.split("..")[0].replace("complement(", "").replace("(", ""))
                    if ira_left_site <= left_site <= ira_right_site:
                        for gene in gene_names:
                            if gene not in ir_gene_list:
                                ir_gene_list.append(gene)
                except ValueError:
                    pass  # Handle cases where location parsing fails
            else:
                judge = False
                joined_parts = temp.replace("join(", "").replace("complement(", "").replace(")", "").split(",")
                for part in joined_parts:
                    try:
                        left_site = int(part.split("..")[0])
                        if ira_left_site <= left_site <= ira_right_site:
                            judge = True
                            break
                    except ValueError:
                        continue  # Skip parts that cannot be parsed
                if judge:
                    for gene in gene_names:
                        if gene not in ir_gene_list:
                            ir_gene_list.append(gene)
    
    return (len(gene_list), len(cds_list), len(trna_list), len(rrna_list), len(ir_gene_list))

def get_organism_and_assembly_name(record):
    organism = record.annotations.get("organism", "Unknown")
    assembly_name = record.annotations.get("comment", "Unknown")
    return organism, assembly_name

def main(): 
    '''
    Main function to process GenBank files and generate a CSV summary.
    '''
    fasta_file_list = get_file_list()
    with open("chloroplast_table.csv", "w") as write_file:
        # Write header including ORGANISM and Assembly Name
        write_file.write(
            "Filename,Total cp genome size (bp),Length of inverted repeat region (bp),"
            "Length of large single copy region (bp),Length of small single copy region (bp),"
            "Total GC content (%),GC content of LSC (%),GC content of IR (%),GC content of SSC (%),"
            "Coding size (bp),Nocoding size (bp),Total number of genes,Number of protein encoding genes,"
            "Number of tRNA genes,Number of rRNA genes,Number of genes duplicated in IR,Organism,Assembly Name\n"
        )
        for each_file in fasta_file_list:
            print(f"{each_file} in process")
            try:
                for record in SeqIO.parse(each_file, "gb"):
                    irb_right_site, irb_left_site, ira_right_site, ira_left_site = find_ir_region(record.seq)
                    gen_num, cds_num, trna_num, rrna_num, ir_gene_num = get_gene_num(record, ira_right_site, ira_left_site)
                    Total_cp_genome_size = len(record.seq)
                    Length_of_inverted_repeat_region = ira_right_site - ira_left_site + 1
                    Length_of_large_single_copy_region = ira_left_site - 1
                    Length_of_small_single_copy_region = irb_left_site - ira_right_site - 1
                    Total_GC_content = GC(record.seq)
                    GC_content_of_LSC = GC(record.seq[:ira_left_site - 1])
                    GC_content_of_IR = GC(record.seq[ira_left_site - 1:irb_left_site])
                    GC_content_of_SSC = GC(record.seq[ira_right_site:irb_left_site - 1])
                    Coding_size = get_coding_size(record)
                    Nocoding_size = len(record.seq) - Coding_size

                    # Get ORGANISM and Assembly Name
                    organism = record.annotations.get("organism", "Unknown")
                    assembly_name = get_assembly_name(each_file)
                    
                    # Write data to CSV
                    write_file.write(
                        f"{each_file},{Total_cp_genome_size},{Length_of_inverted_repeat_region},"
                        f"{Length_of_large_single_copy_region},{Length_of_small_single_copy_region},"
                        f"{Total_GC_content:.2f},{GC_content_of_LSC:.2f},{GC_content_of_IR:.2f},"
                        f"{GC_content_of_SSC:.2f},{Coding_size},{Nocoding_size},"
                        f"{gen_num},{cds_num},{trna_num},{rrna_num},{ir_gene_num},{organism},{assembly_name}\n"
                    )
            except Exception as e:
                print(f"Error processing {each_file}: {e}")

if __name__ == "__main__":
    main()