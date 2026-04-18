import multiprocessing
import os
import shutil
import subprocess
import sys


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_2_coalescent/0_raxml/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/8_2_coalescent/0_raxml/output"
CONDA_ENV_NAME = "trinity_env"

INPUT_EXTENSION = ".fasta"
RAXML_EXECUTABLE_NAME = "raxmlHPC-PTHREADS"

PROCESS_COUNT = max(1, multiprocessing.cpu_count() // 4)
THREADS_PER_PROCESS = 1
BOOTSTRAPS = 1
MODEL = "GTRGAMMA"
SEED = 12345
ADDITIONAL_PARAMETERS = []
DRY_RUN = False
# ============================================================


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


def list_input_files(input_dir):
    return sorted(file_name for file_name in os.listdir(input_dir) if file_name.endswith(INPUT_EXTENSION))


def run_single_raxml(task):
    input_file_name, input_dir, output_dir, raxml_executable = task
    input_path = os.path.join(input_dir, input_file_name)
    prefix = os.path.splitext(input_file_name)[0]
    command = [
        raxml_executable,
        "-T",
        str(THREADS_PER_PROCESS),
        "-n",
        prefix,
        "-s",
        input_path,
        "-m",
        MODEL,
        "-p",
        str(SEED),
        "-f",
        "a",
        "-N",
        str(BOOTSTRAPS),
        "-x",
        str(SEED),
        *ADDITIONAL_PARAMETERS,
    ]
    subprocess.run(
        command,
        cwd=output_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=True,
    )
    return os.path.join(output_dir, f"RAxML_bestTree.{prefix}")


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    os.makedirs(output_dir, exist_ok=True)
    if not os.path.isdir(input_dir):
        print(f"错误: 输入目录不存在: {input_dir}")
        sys.exit(1)
    input_files = list_input_files(input_dir)
    if not input_files:
        print(f"未找到后缀为 {INPUT_EXTENSION} 的文件。")
        sys.exit(0)
    try:
        raxml_executable = resolve_executable(RAXML_EXECUTABLE_NAME)
    except FileNotFoundError as exc:
        print(f"错误: {exc}")
        sys.exit(1)
    if DRY_RUN:
        for input_file in input_files:
            print(f"[DRY RUN] Would build tree for: {input_file}")
        sys.exit(0)
    tasks = [(input_file, input_dir, output_dir, raxml_executable) for input_file in input_files]
    if PROCESS_COUNT == 1:
        outputs = [run_single_raxml(task) for task in tasks]
    else:
        with multiprocessing.Pool(processes=min(PROCESS_COUNT, len(tasks))) as pool:
            outputs = pool.map(run_single_raxml, tasks)
    for output in outputs:
        print(f"建树完成: {output}")


if __name__ == "__main__":
    main()
