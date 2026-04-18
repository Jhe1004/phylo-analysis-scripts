import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import glob

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_1_concatenation/4_iqtree/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_1_concatenation/4_iqtree/output"
CONDA_ENV_NAME = "trinity_env"

INPUT_FILES = []
INPUT_EXTENSIONS = [".fasta", ".fas", ".nex"]
AUTO_CONCAT_SMALL_MATRICES = False

IQTREE_EXECUTABLE_NAME = "iqtree2"
MODEL = "MFP"
THREADS = "AUTO"
BOOTSTRAP_REPLICATES = 1000
ALRT_REPLICATES = 1000
SEED = 42
OVERWRITE_EXISTING = False
PROCESS_COUNT = max(1, multiprocessing.cpu_count() // 2)
EXTRA_ARGS = []
DRY_RUN = False
# ============================================================


PARTITION_FILE = "joint_fragments.partition"
CONCATENATED_ALIGNMENT_FILE = "joint_fragments.fasta"


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


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
    lines = result.stdout.strip().splitlines()
    if not lines:
        return None
    candidate = lines[-1].strip()
    if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return None


def resolve_executable(executable_name):
    local_candidate = os.path.join(SCRIPT_DIR, "dependencies", "bin", executable_name)
    if os.path.isfile(local_candidate) and os.access(local_candidate, os.X_OK):
        return local_candidate
    conda_candidate = find_executable_in_conda_env(executable_name, CONDA_ENV_NAME)
    if conda_candidate:
        return conda_candidate
    path_candidate = shutil.which(executable_name)
    if path_candidate:
        return path_candidate
    raise FileNotFoundError(f"未找到可执行程序: {executable_name}")


def sanitize_taxon_name(name):
    new_name = re.sub(r"[^a-zA-Z0-9_]", "_", name)
    if new_name and new_name[0].isdigit():
        new_name = "T" + new_name
    return new_name[:99]


def discover_input_files(extensions, search_dir):
    discovered = []
    for ext in extensions:
        normalized = ext if ext.startswith(".") else f".{ext}"
        discovered.extend(glob.glob(os.path.join(search_dir, f"*{normalized}")))
    return sorted(set(discovered))


def is_fasta_like_file(filepath):
    return os.path.splitext(filepath)[1].lower() in {".fasta", ".fas", ".fa"}


def create_concatenated_alignment_and_partition(input_files, output_dir):
    fragment_data = []
    taxa_sets = []
    for path in input_files:
        records = {sanitize_taxon_name(r.id): r for r in SeqIO.parse(path, "fasta")}
        if records:
            fragment_data.append((os.path.basename(path), records))
            taxa_sets.append(set(records.keys()))

    if not taxa_sets:
        return None, None

    common_taxa = set.intersection(*taxa_sets)
    if not common_taxa:
        return None, None

    sorted_taxa = sorted(common_taxa)
    concatenated = {taxon: "" for taxon in sorted_taxa}
    partitions = []
    current_pos = 1

    for filename, records in fragment_data:
        frag_name = sanitize_taxon_name(os.path.splitext(filename)[0])
        example_taxon = next(iter(common_taxa))
        frag_len = len(records[example_taxon].seq)
        for taxon in sorted_taxa:
            concatenated[taxon] += str(records[taxon].seq)
        end_pos = current_pos + frag_len - 1
        partitions.append((frag_name, current_pos, end_pos))
        current_pos = end_pos + 1

    concat_path = os.path.join(output_dir, CONCATENATED_ALIGNMENT_FILE)
    partition_path = os.path.join(output_dir, PARTITION_FILE)

    concat_records = []
    for taxon in sorted_taxa:
        concat_records.append(SeqRecord(Seq(concatenated[taxon]), id=taxon, description=""))
    SeqIO.write(concat_records, concat_path, "fasta")

    with open(partition_path, "w", encoding="utf-8") as handle:
        for part_name, start, end in partitions:
            handle.write(f"DNA, {part_name} = {start}-{end}\n")

    return concat_path, partition_path


def build_command(iqtree_executable, input_path, prefix, partition_file=None):
    cmd = [iqtree_executable, "-s", input_path, "-m", str(MODEL), "-T", str(THREADS), "-pre", prefix]
    if partition_file:
        cmd.extend(["-p", partition_file])
    if BOOTSTRAP_REPLICATES:
        cmd.extend(["-B", str(BOOTSTRAP_REPLICATES)])
    if ALRT_REPLICATES:
        cmd.extend(["-alrt", str(ALRT_REPLICATES)])
    if SEED is not None:
        cmd.extend(["-seed", str(SEED)])
    if OVERWRITE_EXISTING:
        cmd.append("-redo")
    cmd.extend([str(x) for x in EXTRA_ARGS])
    return cmd


def run_single_iqtree(task):
    input_path, prefix_name, output_dir, iqtree_executable, partition_file = task
    prefix = os.path.join(output_dir, prefix_name)
    cmd = build_command(iqtree_executable, input_path, prefix, partition_file=partition_file)
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
    return prefix


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    os.makedirs(output_dir, exist_ok=True)
    if not os.path.isdir(input_dir):
        print(f"错误: 输入目录不存在: {input_dir}")
        sys.exit(1)
    input_files = [os.path.join(input_dir, x) for x in INPUT_FILES] if INPUT_FILES else discover_input_files(INPUT_EXTENSIONS, input_dir)
    if not input_files:
        print(f"未找到符合后缀 {INPUT_EXTENSIONS} 的输入文件。")
        sys.exit(0)
    try:
        iqtree_executable = resolve_executable(IQTREE_EXECUTABLE_NAME)
    except FileNotFoundError as exc:
        print(f"错误: {exc}")
        sys.exit(1)

    jobs = []
    fasta_inputs = [path for path in input_files if is_fasta_like_file(path)]
    if AUTO_CONCAT_SMALL_MATRICES and len(fasta_inputs) > 1:
        concat_path, partition_path = create_concatenated_alignment_and_partition(fasta_inputs, output_dir)
        if concat_path and partition_path:
            jobs.append((concat_path, os.path.splitext(os.path.basename(concat_path))[0], output_dir, iqtree_executable, partition_path))

    for input_file in input_files:
        prefix_name = os.path.splitext(os.path.basename(input_file))[0]
        jobs.append((input_file, prefix_name, output_dir, iqtree_executable, None))

    if DRY_RUN:
        for job in jobs:
            print(f"[DRY RUN] Would run IQ-TREE for: {job[0]}")
        sys.exit(0)
    if PROCESS_COUNT == 1:
        outputs = [run_single_iqtree(task) for task in jobs]
    else:
        with multiprocessing.Pool(processes=min(PROCESS_COUNT, len(jobs))) as pool:
            outputs = pool.map(run_single_iqtree, jobs)
    for output in outputs:
        print(f"IQ-TREE 完成: {output}")


if __name__ == "__main__":
    main()
