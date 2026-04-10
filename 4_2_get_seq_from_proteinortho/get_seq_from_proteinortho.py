import os
import sys

import matplotlib
import pandas as pd
import seaborn as sns
from Bio import SeqIO

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
sys.path.append(PROJECT_ROOT)


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
CONDA_ENV_NAME = "trinity_env"

PROTEINORTHO_RESULT_FILE = "proteinortho_result.proteinortho.tsv"
HEATMAP_SPECIES_ORDER_FILE = "name_list.txt"

SUMMARY_FILENAME = "selected_orthologs.tsv"
MAX_ALLOW_MISSING_SPECIES_PROPORTION = 0.5
MAX_ALLOW_LOW_COPY_NUM = 2
MIN_ALG_CONN = 0.1

GENERATE_HEATMAP = False
RANDOM_SEED = 42
HEATMAP_FILENAME = "heatmap.png"

NUCLEOTIDE_OUTPUT_SUFFIX = "_cds.fasta"
PROTEIN_OUTPUT_SUFFIX = "_pep.fasta"
# ============================================================


def resolve_path(*parts):
    return os.path.join(SCRIPT_DIR, *parts)


def validate_paths(input_dir, proteinortho_file):
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not os.path.isfile(proteinortho_file):
        raise FileNotFoundError(f"ProteinOrtho result file not found: {proteinortho_file}")


def load_sequence_files(sequence_dir):
    print("正在将序列文件加载到内存中...")
    file_dict = {}

    for file_name in sorted(os.listdir(sequence_dir)):
        if file_name.endswith(".pep") or file_name.endswith(".cds"):
            file_path = os.path.join(sequence_dir, file_name)
            try:
                file_dict[file_name] = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
            except Exception as exc:
                print(f"无法解析文件 {file_name}: {exc}")

    if not file_dict:
        raise FileNotFoundError(f"在目录 {sequence_dir} 中未找到可解析的 .pep 或 .cds 文件。")

    print("序列文件加载完毕。")
    return file_dict


def get_nucleotide(file_dict, identifier):
    try:
        cds_file = identifier.split("++")[0][:-3] + "cds"
        seq_id = identifier.split("++")[1]
        return str(file_dict[cds_file][seq_id].seq)
    except KeyError:
        return ""


def get_protein(file_dict, identifier):
    try:
        pep_file = identifier.split("++")[0][:-3] + "pep"
        seq_id = identifier.split("++")[1]
        return str(file_dict[pep_file][seq_id].seq)
    except KeyError:
        return ""


def choose_sequence_id(row_value, species_name, file_dict):
    if row_value == "*":
        return None, True

    sequence_ids = row_value.split(",")
    if len(sequence_ids) == 1:
        return sequence_ids[0], True

    if len(sequence_ids) > MAX_ALLOW_LOW_COPY_NUM:
        return None, False

    best_seq_id = ""
    best_length = -1
    for seq_id in sequence_ids:
        full_id = f"{species_name}++{seq_id}"
        sequence = get_nucleotide(file_dict, full_id)
        if len(sequence) > best_length:
            best_length = len(sequence)
            best_seq_id = seq_id
    return best_seq_id, True


def write_ortholog_sequences(selected_row, row_index, columns_values, file_dict, output_dir):
    nucleotide_fasta_path = os.path.join(output_dir, f"ortho{row_index}{NUCLEOTIDE_OUTPUT_SUFFIX}")
    protein_fasta_path = os.path.join(output_dir, f"ortho{row_index}{PROTEIN_OUTPUT_SUFFIX}")

    sequences_to_write = []
    for col_idx in range(3, len(selected_row)):
        species_name = columns_values[col_idx]
        row_value = selected_row.iloc[col_idx]
        chosen_seq_id, is_valid = choose_sequence_id(row_value, species_name, file_dict)
        if not is_valid:
            return False
        if chosen_seq_id:
            sequences_to_write.append((species_name, chosen_seq_id))

    with open(nucleotide_fasta_path, "w", encoding="utf-8") as nucleotide_fasta, open(
        protein_fasta_path, "w", encoding="utf-8"
    ) as protein_fasta:
        for species_name, seq_id in sequences_to_write:
            full_id = f"{species_name}++{seq_id}"
            nucleotide_seq = get_nucleotide(file_dict, full_id)
            protein_seq = get_protein(file_dict, full_id)

            if nucleotide_seq:
                nucleotide_fasta.write(f">{species_name}\n{nucleotide_seq}\n")
            if protein_seq:
                protein_fasta.write(f">{species_name}\n{protein_seq}\n")

    return True


def build_presence_absence_matrix(summary_df):
    species_df = summary_df.iloc[:, 3:].copy()
    return species_df.applymap(lambda value: 0 if value == "*" else 1).T


def get_species_names(df):
    return list(df.columns[3:])


def ensure_species_order_file(species_names):
    order_file = resolve_path(INPUT_DIRECTORY, HEATMAP_SPECIES_ORDER_FILE)
    if not os.path.exists(order_file):
        with open(order_file, "w", encoding="utf-8") as handle:
            for species_name in species_names:
                handle.write(f"{species_name}\n")
        print(f"未检测到 {HEATMAP_SPECIES_ORDER_FILE}，已根据样品自动创建: {order_file}")
    return order_file


def sort_species_for_heatmap(matrix_df):
    order_file = ensure_species_order_file(sorted(matrix_df.index.tolist()))
    with open(order_file, "r", encoding="utf-8") as handle:
        species_order = [line.strip() for line in handle if line.strip()]

    valid_species = [species for species in species_order if species in matrix_df.index]
    remaining_species = [species for species in matrix_df.index if species not in valid_species]
    return matrix_df.loc[valid_species + remaining_species]


def generate_heatmap(summary_df, output_dir):
    print("开始生成热图...")
    matrix_df = build_presence_absence_matrix(summary_df)
    matrix_df = sort_species_for_heatmap(matrix_df)
    matrix_df = matrix_df.sample(frac=1, axis=1, random_state=RANDOM_SEED)

    figure, axis = plt.subplots(figsize=(60, 40))
    sns.heatmap(matrix_df, ax=axis, cmap="GnBu")
    heatmap_path = os.path.join(output_dir, HEATMAP_FILENAME)
    plt.savefig(heatmap_path, dpi=120)
    plt.close()
    print(f"热图已保存: {heatmap_path}")


def main():
    input_dir = resolve_path(INPUT_DIRECTORY)
    output_dir = resolve_path(OUTPUT_DIRECTORY)
    proteinortho_file = os.path.join(input_dir, PROTEINORTHO_RESULT_FILE)
    sequence_dir = input_dir
    os.makedirs(output_dir, exist_ok=True)

    print(f"输入目录: {input_dir}")
    print(f"输出目录: {output_dir}")
    print(f"Conda 环境参数: {CONDA_ENV_NAME}")

    try:
        validate_paths(input_dir, proteinortho_file)
    except FileNotFoundError as exc:
        print(f"错误: {exc}")
        sys.exit(1)

    try:
        df = pd.read_csv(proteinortho_file, sep="\t", header=0)
        ensure_species_order_file(get_species_names(df))
        file_dict = load_sequence_files(sequence_dir)
    except Exception as exc:
        print(f"错误: {exc}")
        sys.exit(1)

    columns_values = df.columns.values
    processed_count = 0
    selected_rows = []

    print("开始处理 ProteinOrtho 结果并提取直系同源序列...")

    for row_index in range(len(df)):
        row = df.iloc[row_index, :]
        missing_ratio = row[3:].tolist().count("*") / len(row[3:])
        if missing_ratio > MAX_ALLOW_MISSING_SPECIES_PROPORTION:
            continue
        if row["Alg.-Conn."] < MIN_ALG_CONN:
            continue

        is_valid = write_ortholog_sequences(row, row_index, columns_values, file_dict, output_dir)
        if not is_valid:
            continue

        processed_count += 1
        selected_rows.append(row)

    if selected_rows:
        selected_df = pd.DataFrame(selected_rows)
        summary_path = os.path.join(output_dir, SUMMARY_FILENAME)
        selected_df.to_csv(summary_path, sep="\t", index=False)
        print(f"已将筛选出的 {len(selected_rows)} 行记录保存到文件: {summary_path}")

        if GENERATE_HEATMAP:
            generate_heatmap(selected_df, output_dir)
    else:
        print("没有找到符合条件的直系同源基因簇。")
        if GENERATE_HEATMAP:
            print("警告：没有找到符合条件的基因簇，无法生成热图。")

    print(f"处理完成！共提取了 {processed_count} 个符合条件的直系同源基因簇。")


if __name__ == "__main__":
    main()
