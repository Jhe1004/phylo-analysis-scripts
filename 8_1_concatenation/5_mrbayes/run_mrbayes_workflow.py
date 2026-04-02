import concurrent.futures
import datetime
import glob
import os
import re
import shutil
import subprocess
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
CONDA_ENV_NAME = "trinity_env"

INPUT_FILES = []
INPUT_EXTENSIONS = [".fasta", ".fas", ".nex"]
PARTITION_FILE = "joint_fragments.nex"
AUTO_CONCAT_SMALL_MATRICES = False

MODEL_NST = "mixed"
MODEL_RATES = "invgamma"
GENERATIONS = 2000000
SAMPLE_FREQ = 1000
CHAINS = 4
RUNS = 2
BURNIN_FRAC = 0.25
THREADS_PER_RUN = 8
MAX_PARALLEL_JOBS = 2

MRBAYES_EXECUTABLE_NAME = "mb"
MPIRUN_EXECUTABLE_NAME = "mpirun"
REPORT_NAME = "mrbayes_summary_report.md"
DRY_RUN = False
# ============================================================


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


MRBAYES_TEMPLATE = """
begin mrbayes;
log start filename = {base_name}_log.txt;

{charsets}
{partition_def}
{set_partition}
set usebeagle=yes beagleresource=0;

lset nst={nst} rates={rates};

{unlinks}

prset applyto=(all) ratepr=variable;

mcmcp ngen={ngen} printfreq={samplefreq} samplefreq={samplefreq} nchains={nchains} nruns={nruns} savebrlens=yes checkpoint=yes checkfreq=1000 filename={base_name}_run;
mcmc;

sumt conformat=Simple contype=Allcompat relburnin=yes burninfrac={burnin};
sump relburnin=yes burninfrac={burnin};
end;
"""


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


def is_nexus_file(filepath):
    return os.path.splitext(filepath)[1].lower() in {".nex", ".nexus"}


def create_nexus_individual(fasta_file, output_dir):
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}.nex")
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        return None
    for record in records:
        record.annotations["molecule_type"] = "DNA"
        record.id = sanitize_taxon_name(record.id)
        record.description = ""
    SeqIO.write(records, output_file, "nexus")
    mb_block = MRBAYES_TEMPLATE.format(
        base_name=os.path.join(output_dir, base_name),
        charsets="",
        partition_def="",
        set_partition="",
        nst=MODEL_NST,
        rates=MODEL_RATES,
        unlinks="unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);",
        ngen=GENERATIONS,
        samplefreq=SAMPLE_FREQ,
        nchains=CHAINS,
        nruns=RUNS,
        burnin=BURNIN_FRAC,
    )
    with open(output_file, "a", encoding="utf-8") as handle:
        handle.write("\n" + "\n".join(line for line in mb_block.split("\n") if line.strip()) + "\n")
    return output_file


def create_nexus_joint(input_files, output_file):
    fragment_data = []
    taxa_sets = []
    for path in input_files:
        records = {sanitize_taxon_name(r.id): r for r in SeqIO.parse(path, "fasta")}
        if records:
            fragment_data.append((os.path.basename(path), records))
            taxa_sets.append(set(records.keys()))
    if not taxa_sets:
        return None
    common_taxa = set.intersection(*taxa_sets)
    if not common_taxa:
        return None
    sorted_taxa = sorted(common_taxa)
    concatenated = {taxon: "" for taxon in sorted_taxa}
    partitions = []
    current_pos = 1
    for filename, records in fragment_data:
        frag_name = os.path.splitext(filename)[0]
        example_taxon = next(iter(common_taxa))
        frag_len = len(records[example_taxon].seq)
        for taxon in sorted_taxa:
            concatenated[taxon] += str(records[taxon].seq)
        end_pos = current_pos + frag_len - 1
        partitions.append((sanitize_taxon_name(frag_name), current_pos, end_pos))
        current_pos = end_pos + 1
    nexus_records = []
    for taxon in sorted_taxa:
        rec = SeqRecord(Seq(concatenated[taxon]), id=taxon, name=taxon, description="")
        rec.annotations["molecule_type"] = "DNA"
        nexus_records.append(rec)
    SeqIO.write(nexus_records, output_file, "nexus")
    charsets = "".join([f"charset {name} = {start}-{end};\n" for name, start, end in partitions])
    partition_def = f"partition by_fragment = {len(partitions)}: {', '.join(name for name, _, _ in partitions)};"
    mb_block = MRBAYES_TEMPLATE.format(
        base_name=os.path.splitext(output_file)[0],
        charsets=charsets,
        partition_def=partition_def,
        set_partition="set partition = by_fragment;",
        nst=MODEL_NST,
        rates=MODEL_RATES,
        unlinks="unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);",
        ngen=GENERATIONS,
        samplefreq=SAMPLE_FREQ,
        nchains=CHAINS,
        nruns=RUNS,
        burnin=BURNIN_FRAC,
    )
    with open(output_file, "a", encoding="utf-8") as handle:
        handle.write("\n" + mb_block + "\n")
    return output_file


def get_nchar(nex_file):
    with open(nex_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            match = re.search(r"dimensions.*nchar\s*=\s*(\d+)", line, re.IGNORECASE)
            if match:
                return int(match.group(1))
    return 1000


def run_mrbayes_task(task):
    nex_file, mpirun_executable, mrbayes_executable = task
    threads = THREADS_PER_RUN
    if get_nchar(nex_file) < 1000 and threads > 2:
        threads = 2
    cmd = [mpirun_executable, "-np", str(threads), mrbayes_executable, nex_file]
    console_log = f"{nex_file}.console.log"
    with open(console_log, "w", encoding="utf-8") as log_f:
        process = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, check=False)
    return nex_file, process.returncode == 0


def parse_mrbayes_stats(base_name):
    metrics = {"lnl": "N/A", "min_ess": "N/A", "max_psrf": "N/A", "status": "Complete"}
    lstat_file = f"{base_name}.lstat"
    if os.path.exists(lstat_file):
        with open(lstat_file, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                parts = line.split()
                if parts and parts[0] == "all":
                    metrics["lnl"] = parts[1]
                    break
    pstat_file = f"{base_name}.pstat"
    if os.path.exists(pstat_file):
        ess_values = []
        psrf_values = []
        with open(pstat_file, "r", encoding="utf-8", errors="ignore") as handle:
            header = None
            for line in handle:
                parts = line.strip().split()
                if not parts:
                    continue
                if "Parameter" in parts and "minESS" in parts:
                    header = parts
                    idx_ess = header.index("minESS")
                    idx_psrf = header.index("PSRF")
                    continue
                if header and len(parts) >= len(header):
                    try:
                        ess_values.append(float(parts[idx_ess]))
                        psrf_values.append(float(parts[idx_psrf]))
                    except ValueError:
                        continue
        if ess_values:
            metrics["min_ess"] = f"{min(ess_values):.2f}"
        if psrf_values:
            metrics["max_psrf"] = f"{max(psrf_values):.4f}"
    return metrics


def generate_report(results, report_path):
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write("# MrBayes Analysis Report\n\n")
        handle.write(f"Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        handle.write("| Dataset | Log Likelihood (Mean) | Min ESS | Max PSRF | Status |\n")
        handle.write("|---------|-----------------------|---------|----------|--------|\n")
        for item in results:
            status = "Success" if item["success"] else "Failed"
            m = item["metrics"]
            handle.write(f"| {item['name']} | {m['lnl']} | {m['min_ess']} | {m['max_psrf']} | {status} |\n")


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    os.makedirs(output_dir, exist_ok=True)
    if not os.path.isdir(input_dir):
        print(f"错误: 输入目录不存在: {input_dir}")
        sys.exit(1)

    input_files = [os.path.join(input_dir, x) for x in INPUT_FILES] if INPUT_FILES else discover_input_files(INPUT_EXTENSIONS, input_dir)
    if not input_files:
        print("未找到可用输入文件。")
        sys.exit(1)

    try:
        mrbayes_executable = resolve_executable(MRBAYES_EXECUTABLE_NAME)
        mpirun_executable = resolve_executable(MPIRUN_EXECUTABLE_NAME)
    except FileNotFoundError as exc:
        print(f"错误: {exc}")
        sys.exit(1)

    jobs = []
    fasta_inputs = [path for path in input_files if is_fasta_like_file(path)]
    if PARTITION_FILE and len(fasta_inputs) > 1 and AUTO_CONCAT_SMALL_MATRICES:
        joint_nex = create_nexus_joint(fasta_inputs, os.path.join(output_dir, PARTITION_FILE))
        if joint_nex:
            jobs.append(joint_nex)
    for input_file in input_files:
        if is_fasta_like_file(input_file):
            nex = create_nexus_individual(input_file, output_dir)
            if nex:
                jobs.append(nex)
        elif is_nexus_file(input_file):
            copied_path = os.path.join(output_dir, os.path.basename(input_file))
            shutil.copy2(input_file, copied_path)
            jobs.append(copied_path)
    jobs = sorted(set(jobs))
    if not jobs:
        print("没有生成可运行的 Nexus 文件。")
        sys.exit(0)

    if DRY_RUN:
        for job in jobs:
            print(f"[DRY RUN] Would run MrBayes for: {job}")
        sys.exit(0)

    tasks = [(job, mpirun_executable, mrbayes_executable) for job in jobs]
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_PARALLEL_JOBS) as executor:
        future_to_job = {executor.submit(run_mrbayes_task, task): task[0] for task in tasks}
        for future in concurrent.futures.as_completed(future_to_job):
            nex = future_to_job[future]
            success = future.result()[1]
            metrics = parse_mrbayes_stats(os.path.splitext(nex)[0] + "_run")
            results.append({"name": os.path.basename(nex), "success": success, "metrics": metrics})

    generate_report(results, os.path.join(output_dir, REPORT_NAME))
    print(f"MrBayes 流程完成，报告已写入: {os.path.join(output_dir, REPORT_NAME)}")


if __name__ == "__main__":
    main()
