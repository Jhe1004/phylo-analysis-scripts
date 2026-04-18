import os
import sys


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_1_concatenation/0_concat/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_1_concatenation/0_concat/output"

INPUT_EXTENSION = ".fas"
OUTPUT_FILENAME = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_1_concatenation/0_concat/output/result.fasta"
DRY_RUN = False
# ============================================================


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def fasta_to_dict(fasta_file):
    result_dict = {}
    current_name = None
    with open(fasta_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                current_name = line.split(" ")[0][1:].strip().replace("/", "_").replace("\\", "_")
                result_dict[current_name] = ""
            else:
                result_dict[current_name] += line.strip()
    any_name = next(iter(result_dict))
    return result_dict, len(result_dict[any_name])


def list_input_files(input_dir):
    return sorted(file_name for file_name in os.listdir(input_dir) if file_name.endswith(INPUT_EXTENSION))


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.isdir(input_dir):
        print(f"错误: 输入目录不存在: {input_dir}")
        sys.exit(1)

    fasta_file_list = list_input_files(input_dir)
    if not fasta_file_list:
        print(f"未找到后缀为 {INPUT_EXTENSION} 的输入文件。")
        sys.exit(0)

    print(f"输入目录: {input_dir}")
    print(f"输出目录: {output_dir}")
    print(f"共找到 {len(fasta_file_list)} 个文件。")

    if DRY_RUN:
        for file_name in fasta_file_list:
            print(f"[DRY RUN] Would concatenate: {file_name}")
        sys.exit(0)

    concat_list = []
    name_list = []
    for file_name in fasta_file_list:
        fasta_dict = fasta_to_dict(os.path.join(input_dir, file_name))
        concat_list.append(fasta_dict)
        for seq_name in fasta_dict[0]:
            if seq_name not in name_list:
                name_list.append(seq_name)

    output_path = os.path.join(output_dir, OUTPUT_FILENAME)
    with open(output_path, "w", encoding="utf-8") as handle:
        for seq_name in name_list:
            handle.write(f">{seq_name}\n")
            for fasta_dict in concat_list:
                if seq_name in fasta_dict[0]:
                    handle.write(fasta_dict[0][seq_name])
                else:
                    handle.write("?" * fasta_dict[1])
            handle.write("\n")

    print(f"拼接完成: {output_path}")


if __name__ == "__main__":
    main()
