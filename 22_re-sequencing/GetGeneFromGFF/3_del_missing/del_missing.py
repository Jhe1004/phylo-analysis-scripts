import os

def filter_fasta_files(folder_path):
    for file in os.listdir(folder_path):
        if file.endswith('.fasta'):
            file_path = os.path.join(folder_path, file)
            filtered_sequences = []
            
            # 读取并过滤序列
            with open(file_path, 'r') as fasta:
                for line in fasta:
                    if line.startswith('>'):
                        seq_name = line.strip()
                    else:
                        seq = line.strip()
                        if seq.count('?') / len(seq) <= 0.8:  # 序列中"?"的比例不超过80%
                            filtered_sequences.append((seq_name, seq))
            
            # 检查过滤后的序列数量
            if len(filtered_sequences) >= 4:
                # 写入过滤后的序列到新文件
                with open(file_path, 'w') as output:
                    for seq_name, seq in filtered_sequences:
                        output.write(f'{seq_name}\n{seq}\n')
            else:
                # 删除序列数量少于4的文件
                os.remove(file_path)

if __name__ == '__main__':
    folder_path = './'
    filter_fasta_files(folder_path)
    print('FASTA files filtered successfully!')
