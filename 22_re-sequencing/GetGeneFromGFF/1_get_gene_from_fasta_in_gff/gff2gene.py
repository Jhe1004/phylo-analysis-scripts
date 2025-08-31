import os

def extract_genes_from_fasta(gff_file, fasta_file, output_file):
    # 读取GFF文件，提取mRNA相关的行，并存储起始和结束位置
    gene_positions = []
    with open(gff_file, 'r') as gff:
        for line in gff:
            if 'mRNA' in line:
                fields = line.strip().split('\t')
                start, end = int(fields[3]), int(fields[4])
                if end - start + 1 <= 10000:  # 只保留长度不超过10000的基因
                    gene_positions.append((start, end))

    # 读取FASTA文件，提取对应位置的序列
    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as output:
        fasta.readline()  # 跳过第一行的标题
        sequence = fasta.read().replace('\n', '')
        for i, (start, end) in enumerate(gene_positions, start=1):
            gene_sequence = sequence[start-1:end]
            output.write(f'>gene_{i}\n{gene_sequence}\n')

def process_all_fasta_files(folder_path, gff_file):
    for file in os.listdir(folder_path):
        if file.endswith('.fasta'):
            fasta_file = os.path.join(folder_path, file)
            output_file = os.path.join(folder_path, f'extracted_{file}')
            extract_genes_from_fasta(gff_file, fasta_file, output_file)
            print(f'Processed {fasta_file}')

if __name__ == '__main__':
    folder_path = './'
    gff_file = 'ref.gff'
    process_all_fasta_files(folder_path, gff_file)
    print('All files processed successfully!')
