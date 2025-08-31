import os
from collections import defaultdict

def extract_sequences_to_files(folder_path, output_folder):
    sequences = defaultdict(list)
    
    # 遍历文件夹中的所有FASTA文件
    for file in os.listdir(folder_path):
        if file.endswith('.fasta'):
            with open(os.path.join(folder_path, file), 'r') as fasta:
                for line in fasta:
                    if line.startswith('>'):
                        seq_name = line.strip().replace('>', '')
                    else:
                        sequences[seq_name].append(str(file) + "~" + line.strip())
    
    # 为每个序列名称创建一个新的FASTA文件
    for seq_name, seq_list in sequences.items():
        with open(os.path.join(output_folder, f'{seq_name}.fasta'), 'w') as output:
            for seq in seq_list:
                output.write(f'>{seq.split("~")[0]}\n{seq.split("~")[1]}\n')

if __name__ == '__main__':
    folder_path = './'
    output_folder = './'
    extract_sequences_to_files(folder_path, output_folder)
    print('Sequences extracted to individual files successfully!')