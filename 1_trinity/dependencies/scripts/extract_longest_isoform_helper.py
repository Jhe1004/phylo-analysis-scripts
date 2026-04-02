import glob
import os
import re
import sys

try:
    from Bio import SeqIO
except ImportError:
    print("错误: 缺少 Biopython 库。请在当前 Python 环境中安装 biopython。")
    sys.exit(1)


def get_gene_id(seq_id, regex_pattern):
    match = regex_pattern.search(seq_id)
    if match:
        return match.group(1)

    parts = seq_id.rsplit("_", 1)
    if len(parts) > 1:
        return parts[0]
    return None


def process_single_file(filepath, output_path, regex_pattern, logger):
    gene_dict = {}
    total_count = 0

    logger.info(f"Processing FASTA file: {filepath}")

    for record in SeqIO.parse(filepath, "fasta"):
        total_count += 1
        gene_id = get_gene_id(record.id, regex_pattern)
        if not gene_id:
            logger.warning(f"Cannot parse gene ID, skipped: {record.id}")
            continue

        current_length = len(record.seq)
        if gene_id not in gene_dict or current_length > len(gene_dict[gene_id].seq):
            gene_dict[gene_id] = record

    longest_transcripts = sorted(gene_dict.values(), key=lambda item: len(item.seq), reverse=True)

    if not longest_transcripts:
        logger.warning(f"No valid sequence found in file: {filepath}")
        return

    SeqIO.write(longest_transcripts, output_path, "fasta")
    logger.info(
        f"Saved longest isoform file: {output_path} "
        f"(input sequences: {total_count}, kept genes: {len(longest_transcripts)})"
    )


def run_extract_longest_isoform(
    input_directory,
    output_directory,
    input_file_pattern,
    output_extension,
    gene_id_regex,
    logger,
):
    if not os.path.isdir(input_directory):
        raise FileNotFoundError(f"Input directory not found: {input_directory}")

    regex_pattern = re.compile(gene_id_regex)
    search_pattern = os.path.join(input_directory, input_file_pattern)
    input_files = sorted(glob.glob(search_pattern))
    os.makedirs(output_directory, exist_ok=True)

    if not input_files:
        logger.warning(f"No files matched pattern: {search_pattern}")
        return

    logger.info(f"Found {len(input_files)} Trinity FASTA files for longest isoform extraction")

    for input_file in input_files:
        sample_name = os.path.basename(os.path.dirname(input_file))
        output_file = f"{sample_name}{output_extension}"
        output_path = os.path.join(output_directory, output_file)
        process_single_file(input_file, output_path, regex_pattern, logger)
