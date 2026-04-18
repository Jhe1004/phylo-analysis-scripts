import multiprocessing
import os
import shutil
import sys
from functools import partial

import pandas as pd
from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
sys.path.append(PROJECT_ROOT)

try:
    from lib import utils
except ImportError:
    print("Error: Could not import 'lib'. Please ensure the project structure is intact.")
    sys.exit(1)


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/6_2_del_indel/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/6_2_del_indel/output"
CONDA_ENV_NAME = "trinity_env"

INPUT_EXTENSION = ".fasta"
OUTPUT_EXTENSION = ".fas"

ALIGNMENT_LEN = 50000
MISSING_SITE_PROPORTION = 0.2
PROCESS_COUNT = max(1, multiprocessing.cpu_count() // 2)

NORMAL_SITE_LIST = ["A", "a", "T", "t", "C", "c", "G", "g"]
TEMP_DIRECTORY_NAME = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/6_2_del_indel/output/temp"
SPLIT_TAG = ".split."
DRY_RUN = False
# ============================================================


def validate_config(input_dir):
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if MISSING_SITE_PROPORTION == 0:
        raise ValueError("MISSING_SITE_PROPORTION cannot be 0.")


def list_input_files(input_dir):
    return sorted(file_name for file_name in os.listdir(input_dir) if file_name.endswith(INPUT_EXTENSION))


def sanitize_and_pad_fasta(input_path, output_path):
    valid_records = []
    max_len = 0

    for record in SeqIO.parse(input_path, "fasta"):
        seq_str = str(record.seq)
        if not all(char in "-?Nn" for char in seq_str):
            valid_records.append(record)
            max_len = max(max_len, len(seq_str))

    if not valid_records:
        return False, 0

    with open(output_path, "w", encoding="utf-8") as handle:
        for record in valid_records:
            seq_str = str(record.seq)
            if len(seq_str) < max_len:
                seq_str += "-" * (max_len - len(seq_str))
            handle.write(f">{record.id}\n{seq_str}\n")

    return True, max_len


def split_fasta_file(fasta_path, alignment_len):
    generated_files = []
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return generated_files

    sequence_length = len(records[0].seq)
    left = 0
    while left < sequence_length:
        right = min(left + alignment_len, sequence_length)
        split_path = f"{fasta_path[:-3]}{SPLIT_TAG}{left}.fa"
        with open(split_path, "w", encoding="utf-8") as handle:
            for record in records:
                subseq = str(record.seq[left:right])
                handle.write(f">{record.id}\n{subseq}\n")
        generated_files.append(split_path)
        left += alignment_len
    os.remove(fasta_path)
    return generated_files


def preprocess_single_file(input_file_name, input_dir, temp_dir):
    input_path = os.path.join(input_dir, input_file_name)
    base_name = input_file_name[: -len(INPUT_EXTENSION)]
    normalized_fa = os.path.join(temp_dir, f"{base_name}.fa")

    is_valid, max_len = sanitize_and_pad_fasta(input_path, normalized_fa)
    if not is_valid:
        return {"input_file": input_file_name, "generated_files": [], "was_split": False, "is_valid": False}

    if max_len >= ALIGNMENT_LEN:
        generated_files = split_fasta_file(normalized_fa, ALIGNMENT_LEN)
        return {"input_file": input_file_name, "generated_files": generated_files, "was_split": True, "is_valid": True}

    return {"input_file": input_file_name, "generated_files": [normalized_fa], "was_split": False, "is_valid": True}


def get_seq_name_list(fasta_path):
    seq_name_list = []
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                seq_name_list.append(line)
    return seq_name_list


def calculate_filtered_sequences(fasta_path):
    seq_array_temp_list = []
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                seq_array_temp_list.append(list(line.strip()))

    seq_array = pd.DataFrame(seq_array_temp_list)
    if seq_array.empty:
        return []

    row_count = seq_array.shape[0]
    kept_columns = []
    for column_index in range(seq_array.shape[1]):
        gap_num = 0
        for site in list(seq_array[column_index]):
            if site not in NORMAL_SITE_LIST:
                gap_num += 1
        if gap_num / row_count < MISSING_SITE_PROPORTION:
            kept_columns.append(list(seq_array[column_index]))

    temp_seq_array = pd.DataFrame(kept_columns)
    result_list = []
    for row_index in range(row_count):
        try:
            result_list.append("".join(list(temp_seq_array[row_index].values)))
        except Exception:
            return []
    return result_list


def write_filtered_file(output_path, seq_name_list, seq_list):
    with open(output_path, "w", encoding="utf-8") as handle:
        for index in range(len(seq_name_list)):
            judge1 = seq_list[index].count("-") == len(seq_list[index])
            judge2 = seq_list[index].count("?") == len(seq_list[index])
            if judge1 or judge2:
                continue
            handle.write(seq_name_list[index])
            handle.write(seq_list[index] + "\n")


def filter_single_alignment(fasta_path):
    seq_name_list = get_seq_name_list(fasta_path)
    seq_list = calculate_filtered_sequences(fasta_path)
    if not seq_list:
        return None

    output_path = fasta_path + "s"
    write_filtered_file(output_path, seq_name_list, seq_list)
    os.remove(fasta_path)
    return output_path


def fasta_to_dict(fasta_file):
    result_dict = {}
    current_seq_name = None
    with open(fasta_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                current_seq_name = line.split(" ")[0][1:].strip().replace("/", "_").replace("\\", "_")
                result_dict[current_seq_name] = ""
            else:
                result_dict[current_seq_name] += line.strip()
    any_seq_name = next(iter(result_dict))
    return result_dict, len(result_dict[any_seq_name])


def concat_split_gene(gene_name, split_output_files, output_dir):
    fasta_file_list = sorted(path for path in split_output_files if os.path.basename(path).startswith(f"{gene_name}{SPLIT_TAG}"))
    concat_list = []
    name_list = []

    for fasta_file in fasta_file_list:
        fasta_dict = fasta_to_dict(fasta_file)
        concat_list.append(fasta_dict)
        for each_name in fasta_dict[0]:
            if each_name not in name_list:
                name_list.append(each_name)

    output_path = os.path.join(output_dir, f"{gene_name}{OUTPUT_EXTENSION}")
    with open(output_path, "w", encoding="utf-8") as handle:
        for each_name in name_list:
            handle.write(f">{each_name}\n")
            for each_dict in concat_list:
                if each_name in each_dict[0]:
                    handle.write(each_dict[0][each_name])
                else:
                    handle.write("?" * each_dict[1])
            handle.write("\n")

    for fasta_file in fasta_file_list:
        os.remove(fasta_file)

    return output_path


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    temp_dir = os.path.join(output_dir, TEMP_DIRECTORY_NAME)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    logger = utils.setup_logger("DelIndel", os.path.join(output_dir, "del_indel.log"))

    logger.info("Starting delete-indel workflow")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Default conda environment: {CONDA_ENV_NAME}")
    logger.info(f"Using multiprocessing: process_count={PROCESS_COUNT}")

    try:
        validate_config(input_dir)
    except (FileNotFoundError, ValueError) as exc:
        logger.error(str(exc))
        sys.exit(1)

    input_files = list_input_files(input_dir)
    if not input_files:
        logger.warning(f"No input files found with extension {INPUT_EXTENSION}")
        sys.exit(0)

    if DRY_RUN:
        for input_file in input_files:
            logger.info(f"[DRY RUN] Would process: {input_file}")
        sys.exit(0)

    preprocess_worker = partial(preprocess_single_file, input_dir=input_dir, temp_dir=temp_dir)
    if PROCESS_COUNT == 1:
        preprocess_results = [preprocess_worker(input_file) for input_file in input_files]
    else:
        with multiprocessing.Pool(processes=min(PROCESS_COUNT, len(input_files))) as pool:
            preprocess_results = pool.map(preprocess_worker, input_files)

    fa_files = []
    split_gene_map = {}
    for result in preprocess_results:
        if not result["is_valid"]:
            logger.warning(f"No valid samples remained after preprocessing: {result['input_file']}")
            continue
        fa_files.extend(result["generated_files"])
        if result["was_split"]:
            gene_name = result["input_file"][: -len(INPUT_EXTENSION)]
            split_gene_map[gene_name] = []

    if not fa_files:
        logger.warning("No valid intermediate alignment files generated.")
        sys.exit(0)

    if PROCESS_COUNT == 1:
        filtered_outputs = [filter_single_alignment(fa_file) for fa_file in fa_files]
    else:
        with multiprocessing.Pool(processes=min(PROCESS_COUNT, len(fa_files))) as pool:
            filtered_outputs = pool.map(filter_single_alignment, fa_files)

    filtered_outputs = [path for path in filtered_outputs if path]
    if not filtered_outputs:
        logger.warning("No filtered alignment files were produced.")
        sys.exit(0)

    for output_path in filtered_outputs:
        file_name = os.path.basename(output_path)
        if SPLIT_TAG in file_name:
            gene_name = file_name.split(SPLIT_TAG)[0]
            split_gene_map.setdefault(gene_name, []).append(output_path)
        else:
            final_path = os.path.join(output_dir, file_name)
            shutil.move(output_path, final_path)
            logger.info(f"Finished: {final_path}")

    if split_gene_map:
        concat_worker = partial(concat_split_gene, split_output_files=filtered_outputs, output_dir=output_dir)
        gene_names = sorted(split_gene_map.keys())
        if PROCESS_COUNT == 1:
            concat_outputs = [concat_worker(gene_name) for gene_name in gene_names]
        else:
            with multiprocessing.Pool(processes=min(PROCESS_COUNT, len(gene_names))) as pool:
                concat_outputs = pool.map(concat_worker, gene_names)
        for concat_output in concat_outputs:
            logger.info(f"Concatenated split alignment: {concat_output}")

    shutil.rmtree(temp_dir, ignore_errors=True)
    logger.info("Delete-indel workflow finished")


if __name__ == "__main__":
    main()
