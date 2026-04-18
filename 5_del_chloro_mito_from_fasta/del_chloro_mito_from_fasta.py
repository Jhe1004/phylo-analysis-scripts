import multiprocessing
import os
import shutil
import subprocess
import sys
from functools import partial

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
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/5_del_chloro_mito_from_fasta/output"
CONDA_ENV_NAME = "trinity_env"

REFERENCE_GENBANK_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/5_del_chloro_mito_from_fasta/reference_gb"
SEQUENCE_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/5_del_chloro_mito_from_fasta/sequences"

PEP_SUFFIX = "pep.fasta"
CDS_SUFFIX = "cds.fasta"

MAKEBLASTDB_EXECUTABLE_NAME = "makeblastdb"
BLASTP_EXECUTABLE_NAME = "blastp"

PROCESS_COUNT = min(multiprocessing.cpu_count(), 64)
BLAST_THREADS_PER_PROCESS = 1
BLAST_TASK = "blastp-short"
EVALUE = "0.00001"
CONTAMINATION_RATIO_THRESHOLD = 0.1

DRY_RUN = False
# ============================================================


def validate_paths(reference_dir, sequence_dir):
    if not os.path.isdir(reference_dir):
        raise FileNotFoundError(f"Reference genbank directory not found: {reference_dir}")
    if not os.path.isdir(sequence_dir):
        raise FileNotFoundError(f"Sequence directory not found: {sequence_dir}")


def find_executable_in_conda_env(executable_name, env_name):
    try:
        result = subprocess.run(
            ["conda", "run", "-n", env_name, "which", executable_name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=True,
        )
    except Exception:
        return None

    resolved_path = result.stdout.strip().splitlines()
    if not resolved_path:
        return None

    candidate = resolved_path[-1].strip()
    if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return None


def resolve_executable(executable_name, logger):
    local_candidate = os.path.join(SCRIPT_DIR, "dependencies", "bin", executable_name)
    if os.path.isfile(local_candidate) and os.access(local_candidate, os.X_OK):
        logger.info(f"Using local dependency executable: {local_candidate}")
        return local_candidate

    conda_candidate = find_executable_in_conda_env(executable_name, CONDA_ENV_NAME)
    if conda_candidate:
        logger.info(f"Using executable from conda env '{CONDA_ENV_NAME}': {conda_candidate}")
        return conda_candidate

    path_candidate = shutil.which(executable_name)
    if path_candidate:
        logger.info(f"Using executable from PATH: {path_candidate}")
        return path_candidate

    raise FileNotFoundError(
        f"Required executable not found: {executable_name}. "
        f"Searched in dependencies/bin, conda env '{CONDA_ENV_NAME}', and PATH."
    )


def list_reference_gb_files(reference_dir):
    return sorted(
        os.path.join(reference_dir, file_name)
        for file_name in os.listdir(reference_dir)
        if file_name.endswith(".gb")
    )


def list_pep_files(sequence_dir):
    return sorted(
        file_name for file_name in os.listdir(sequence_dir) if file_name.endswith(PEP_SUFFIX)
    )


def extract_reference_proteins(reference_gb_files, reference_protein_path):
    seen_genes = set()

    with open(reference_protein_path, "w", encoding="utf-8") as handle:
        for gb_file in reference_gb_files:
            for record in SeqIO.parse(gb_file, "genbank"):
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    gene_values = feature.qualifiers.get("gene")
                    translation_values = feature.qualifiers.get("translation")
                    if not gene_values or not translation_values:
                        continue
                    gene_name = gene_values[0]
                    protein_seq = translation_values[0]
                    if gene_name in seen_genes:
                        continue
                    seen_genes.add(gene_name)
                    handle.write(f">{gene_name}\n{protein_seq}\n")


def count_fasta_records(fasta_path):
    count = 0
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def build_blast_db(makeblastdb_executable, reference_protein_path, db_prefix, logger):
    cmd = [
        makeblastdb_executable,
        "-in",
        reference_protein_path,
        "-input_type",
        "fasta",
        "-dbtype",
        "prot",
        "-out",
        db_prefix,
        "-title",
        os.path.basename(db_prefix),
    ]
    utils.run_command(cmd, logger, cwd=os.path.dirname(reference_protein_path), dry_run=DRY_RUN)


def analyse_single_pep_file(
    pep_file_name,
    sequence_dir,
    db_prefix,
    blastp_executable,
    temp_dir,
):
    pep_path = os.path.join(sequence_dir, pep_file_name)
    total_sequences = count_fasta_records(pep_path)

    output_file = os.path.join(temp_dir, f"{pep_file_name}.blast.tsv")
    cmd = [
        blastp_executable,
        "-task",
        BLAST_TASK,
        "-query",
        pep_path,
        "-db",
        db_prefix,
        "-out",
        output_file,
        "-evalue",
        EVALUE,
        "-outfmt",
        "6 qseqid",
        "-num_threads",
        str(BLAST_THREADS_PER_PROCESS),
    ]

    subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=True,
    )

    hit_query_ids = set()
    if os.path.exists(output_file):
        with open(output_file, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if line:
                    hit_query_ids.add(line.split("\t")[0])

    hit_count = len(hit_query_ids)
    contamination_ratio = (hit_count / total_sequences) if total_sequences else 0.0
    is_contaminated = contamination_ratio >= CONTAMINATION_RATIO_THRESHOLD

    return {
        "pep_file": pep_file_name,
        "total_sequences": total_sequences,
        "hit_count": hit_count,
        "contamination_ratio": contamination_ratio,
        "is_contaminated": is_contaminated,
    }


def copy_result_files(result, sequence_dir, kept_dir, removed_dir):
    pep_file = result["pep_file"]
    cds_file = pep_file.replace("pep", "cds")

    source_pep = os.path.join(sequence_dir, pep_file)
    source_cds = os.path.join(sequence_dir, cds_file)

    target_root = removed_dir if result["is_contaminated"] else kept_dir
    shutil.copyfile(source_pep, os.path.join(target_root, pep_file))
    if os.path.exists(source_cds):
        shutil.copyfile(source_cds, os.path.join(target_root, cds_file))


def write_summary(results, output_dir):
    summary_path = os.path.join(output_dir, "summary.tsv")
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write("pep_file\ttotal_sequences\thit_count\tcontamination_ratio\tis_contaminated\n")
        for result in results:
            handle.write(
                f"{result['pep_file']}\t{result['total_sequences']}\t{result['hit_count']}\t"
                f"{result['contamination_ratio']:.6f}\t{result['is_contaminated']}\n"
            )

    removed_list_path = os.path.join(output_dir, "removed_files.txt")
    with open(removed_list_path, "w", encoding="utf-8") as handle:
        for result in results:
            if result["is_contaminated"]:
                pep_file = result["pep_file"]
                cds_file = pep_file.replace("pep", "cds")
                handle.write(f"{pep_file}\n{cds_file}\n")


def main():
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    reference_dir = os.path.join(SCRIPT_DIR, REFERENCE_GENBANK_DIRECTORY)
    sequence_dir = os.path.join(SCRIPT_DIR, SEQUENCE_DIRECTORY)
    db_dir = os.path.join(output_dir, "reference_db")
    temp_dir = os.path.join(output_dir, "temp")
    kept_dir = os.path.join(output_dir, "kept_sequences")
    removed_dir = os.path.join(output_dir, "removed_sequences")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(kept_dir, exist_ok=True)
    os.makedirs(removed_dir, exist_ok=True)
    logger = utils.setup_logger("DelChloroMito", os.path.join(output_dir, "del_chloro_mito.log"))

    logger.info("Starting chloroplast/mitochondria filtering workflow")
    logger.info(f"Reference genbank directory: {reference_dir}")
    logger.info(f"Sequence directory: {sequence_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Default conda environment: {CONDA_ENV_NAME}")

    try:
        validate_paths(reference_dir, sequence_dir)
        makeblastdb_executable = resolve_executable(MAKEBLASTDB_EXECUTABLE_NAME, logger)
        blastp_executable = resolve_executable(BLASTP_EXECUTABLE_NAME, logger)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        sys.exit(1)

    reference_gb_files = list_reference_gb_files(reference_dir)
    if not reference_gb_files:
        logger.error("No reference .gb files found.")
        sys.exit(1)

    pep_files = list_pep_files(sequence_dir)
    if not pep_files:
        logger.warning(f"No peptide files found with suffix: {PEP_SUFFIX}")
        sys.exit(0)

    reference_protein_path = os.path.join(db_dir, "ref.protein")
    db_prefix = os.path.join(db_dir, "ref_db")
    extract_reference_proteins(reference_gb_files, reference_protein_path)
    build_blast_db(makeblastdb_executable, reference_protein_path, db_prefix, logger)

    logger.info(f"Found {len(pep_files)} peptide files to process")
    logger.info(
        f"Using multiprocessing: process_count={PROCESS_COUNT}, "
        f"blast_threads_per_process={BLAST_THREADS_PER_PROCESS}"
    )

    worker = partial(
        analyse_single_pep_file,
        sequence_dir=sequence_dir,
        db_prefix=db_prefix,
        blastp_executable=blastp_executable,
        temp_dir=temp_dir,
    )

    try:
        if DRY_RUN:
            results = []
            for pep_file_name in pep_files:
                logger.info(f"[DRY RUN] Would analyse: {pep_file_name}")
                results.append(
                    {
                        "pep_file": pep_file_name,
                        "total_sequences": 0,
                        "hit_count": 0,
                        "contamination_ratio": 0.0,
                        "is_contaminated": False,
                    }
                )
        elif PROCESS_COUNT == 1:
            results = [worker(pep_file_name) for pep_file_name in pep_files]
        else:
            with multiprocessing.Pool(processes=min(PROCESS_COUNT, len(pep_files))) as pool:
                results = pool.map(worker, pep_files)
    except Exception as exc:
        logger.error(f"Parallel blastp analysis failed: {exc}")
        sys.exit(1)

    for result in results:
        copy_result_files(result, sequence_dir, kept_dir, removed_dir)
        logger.info(
            f"{result['pep_file']}: total={result['total_sequences']}, hits={result['hit_count']}, "
            f"ratio={result['contamination_ratio']:.4f}, removed={result['is_contaminated']}"
        )

    write_summary(results, output_dir)
    logger.info("Filtering workflow finished")


if __name__ == "__main__":
    main()
