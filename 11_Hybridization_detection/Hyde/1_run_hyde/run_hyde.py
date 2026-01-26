# -*- coding: utf-8 -*-
import argparse
import sys
import os

def format_fasta_to_single_line(original_fasta_path):
    """
    读取一个FASTA文件（可以是多行换行的格式），
    并将其转换为每个序列只占一行的标准格式。
    返回新生成的文件名。
    """
    base_name = os.path.splitext(original_fasta_path)[0]
    linear_fasta_path = f"{base_name}.linear.fasta"
    sequences = {}
    current_seq_id = None

    with open(original_fasta_path, "r") as f_in:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_seq_id = line[1:]
                sequences[current_seq_id] = ""
            elif current_seq_id:
                sequences[current_seq_id] += line
    
    with open(linear_fasta_path, "w") as f_out:
        for seq_id, sequence in sequences.items():
            f_out.write(f">{seq_id}\n")
            f_out.write(f"{sequence}\n")
    
    return linear_fasta_path

def fasta2phy(fasta_file):
    """
    将单行格式的fasta序列转化为phy文件，并返回物种数目、位点数目和phy文件名。
    """
    species_num = 0
    species_len = 0
    with open(fasta_file, "r") as read_file:
        for each_line in read_file:
            if len(each_line) > 1:
                if each_line.startswith(">"):
                    species_num += 1
                else:
                    species_len = len(each_line.strip())
    
    phy_file = fasta_file.replace(".fasta", "") + ".phy"
    with open(phy_file, "w") as write_file:
        with open(fasta_file, "r") as read_file:
            write_file.write(f" {species_num} {species_len}\n")
            header = ''
            sequence_started = False
            for each_line in read_file:
                line = each_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if sequence_started:
                        write_file.write("\n")
                    header = line.replace(">", "").strip()
                    write_file.write(f"{header} ")
                    sequence_started = True
                else:
                    write_file.write(line)
            write_file.write("\n")
            
    return str(species_num), str(species_len), phy_file
    
def make_mapfile(fasta_file, outgroup_name, map_filename):
    """
    根据FASTA文件和指定的外群名称，生成指定名称的map文件。
    """
    with open(map_filename, "w") as write_file:
        with open(fasta_file, "r") as read_file:
            for each_line in read_file:
                if each_line.startswith(">"):
                    species_name = each_line.replace(">", "").strip()
                    if species_name == outgroup_name:
                        write_file.write(f"{species_name}\tout\n")
                    else:
                        write_file.write(f"{species_name}\t{species_name}\n")

def main():
    """
    主执行函数
    """
    parser = argparse.ArgumentParser(
        description="Run HyDe in batches. Automatically finds all .fasta files, uses max CPU threads, processes them, and cleans up.",
        add_help=True
    )
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-o', '--outgroup', action="store", metavar='\b', type=str, required=True, 
                          help="Name of the outgroup species (must match a FASTA header exactly)")

    args = parser.parse_args()
    outgroup_name = args.outgroup

    # --- 核心修改：自动检测CPU核心数 ---
    try:
        max_threads = os.cpu_count()
        if max_threads is None:
            print("警告：无法自动检测CPU核心数，将使用默认值 4。")
            max_threads = 4
        else:
            print(f"信息：已自动检测到本机最大线程数为 {max_threads}，将用于 -j 参数。")
    except NotImplementedError:
        print("警告：此平台不支持检测CPU核心数，将使用默认值 4。")
        max_threads = 4
    
    fasta_files_to_process = [f for f in os.listdir('.') if f.endswith('.fasta')]

    if not fasta_files_to_process:
        print("错误：在当前目录下没有找到任何 .fasta 文件。")
        sys.exit(1)

    print(f"检测到 {len(fasta_files_to_process)} 个 FASTA 文件，将开始批量处理...")
    
    for fasta_file in fasta_files_to_process:
        print(f"\n{'='*20}\n[+] 正在处理文件: {fasta_file}\n{'='*20}")
        
        base_name = os.path.splitext(fasta_file)[0]
        linear_fasta_file = f"{base_name}.linear.fasta"
        infile_phy = f"{base_name}.linear.phy"
        map_filename = f"{base_name}_map.txt"
        
        try:
            # 1. 格式化FASTA
            print(f"  (1/5) 格式化FASTA为单行序列...")
            format_fasta_to_single_line(fasta_file)

            # 2. 转换为PHY格式
            print(f"  (2/5) 转换格式为PHYLIP...")
            species_num, species_len, _ = fasta2phy(linear_fasta_file)

            # 3. 创建唯一的Map文件
            print(f"  (3/5) 创建独立的map文件...")
            make_mapfile(linear_fasta_file, outgroup_name, map_filename)

            # 4. 运行HyDe (使用检测到的线程数)
            print(f"  (4/5) 准备并运行HyDe...")
            hyde_command = (
                f"python hyde.py -i {infile_phy} -j {max_threads} --ignore_amb_sites " # 核心修改点
                f"-m {map_filename} -o out -n {species_num} -t {species_num} -s {species_len} "
                f"--prefix {base_name}"
            )
            print(f"      -> 即将执行命令: {hyde_command}")
            os.system(hyde_command)

            # 5. 清理中间文件
            print(f"  (5/5) 清理中间文件...")
            try:
                os.remove(linear_fasta_file)
                os.remove(infile_phy)
                os.remove(map_filename)
            except OSError as e:
                print(f"      -> 清理文件时出错: {e}")

            print(f"\n[✓] 文件 {fasta_file} 处理完成.")

        except Exception as e:
            print(f"\n[!] 处理文件 {fasta_file} 时发生严重错误: {e}")
            continue

    print(f"\n{'='*20}\n所有任务已完成！\n{'='*20}")

if __name__ == "__main__":
    main()