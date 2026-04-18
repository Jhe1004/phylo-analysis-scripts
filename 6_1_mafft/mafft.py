import multiprocessing
import os
import shutil
import subprocess
import sys
from functools import partial

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
sys.path.append(PROJECT_ROOT)

try:
    from lib import utils
except ImportError:
    print("Error: Could not import 'lib'. Please ensure the project structure is intact.")
    sys.exit(1)


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/6_1_mafft/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/6_1_mafft/output"
CONDA_ENV_NAME = "trinity_env"

INPUT_EXTENSION = ".fasta"
OUTPUT_SUFFIX = "_maffted.fasta"

MAFFT_EXECUTABLE_NAME = "mafft"
PROCESS_COUNT = max(1, multiprocessing.cpu_count() // 2)
THREADS_PER_PROCESS = 4

ADJUST_DIRECTION = True
AUTO_MODE = True
REMOVE_REVERSED_TAG = True

ADDITIONAL_PARAMETERS = []
DRY_RUN = False
# ============================================================


def validate_config(input_dir):
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")


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


def list_input_files(input_dir):
    return sorted(file_name for file_name in os.listdir(input_dir) if file_name.endswith(INPUT_EXTENSION))


def clean_sequence_headers(fasta_text):
    if not REMOVE_REVERSED_TAG:
        return fasta_text
    return fasta_text.replace("_R_", "")


def run_single_mafft(input_file_name, input_dir, output_dir, mafft_executable):
    input_path = os.path.join(input_dir, input_file_name)
    base_name = input_file_name[: -len(INPUT_EXTENSION)] if input_file_name.endswith(INPUT_EXTENSION) else input_file_name
    output_path = os.path.join(output_dir, f"{base_name}{OUTPUT_SUFFIX}")

    cmd = [
        mafft_executable,
        "--thread",
        str(THREADS_PER_PROCESS),
        *ADDITIONAL_PARAMETERS,
    ]

    if ADJUST_DIRECTION:
        cmd.append("--adjustdirection")
    if AUTO_MODE:
        cmd.append("--auto")

    cmd.append(input_path)

    process = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=True,
    )

    output_text = clean_sequence_headers(process.stdout)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(output_text)

    return {
        "input_file": input_file_name,
        "output_file": output_path,
        "stderr": process.stderr.strip(),
    }


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    os.makedirs(output_dir, exist_ok=True)
    logger = utils.setup_logger("MAFFT", os.path.join(output_dir, "mafft.log"))

    logger.info("Starting MAFFT workflow")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Default conda environment: {CONDA_ENV_NAME}")
    logger.info(
        f"Using multiprocessing: process_count={PROCESS_COUNT}, "
        f"threads_per_process={THREADS_PER_PROCESS}"
    )

    try:
        validate_config(input_dir)
        mafft_executable = resolve_executable(MAFFT_EXECUTABLE_NAME, logger)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        sys.exit(1)

    input_files = list_input_files(input_dir)
    if not input_files:
        logger.warning(f"No input files found with extension {INPUT_EXTENSION}")
        sys.exit(0)

    logger.info(f"Found {len(input_files)} files to process: {input_files}")

    worker = partial(
        run_single_mafft,
        input_dir=input_dir,
        output_dir=output_dir,
        mafft_executable=mafft_executable,
    )

    try:
        if DRY_RUN:
            for input_file_name in input_files:
                logger.info(f"[DRY RUN] Would process: {input_file_name}")
        elif PROCESS_COUNT == 1:
            results = [worker(input_file_name) for input_file_name in input_files]
            for result in results:
                logger.info(f"Finished: {result['input_file']} -> {result['output_file']}")
                if result["stderr"]:
                    logger.info(f"STDERR:\n{result['stderr']}")
        else:
            with multiprocessing.Pool(processes=min(PROCESS_COUNT, len(input_files))) as pool:
                results = pool.map(worker, input_files)
            for result in results:
                logger.info(f"Finished: {result['input_file']} -> {result['output_file']}")
                if result["stderr"]:
                    logger.info(f"STDERR:\n{result['stderr']}")
    except Exception as exc:
        logger.error(f"MAFFT workflow failed: {exc}")
        sys.exit(1)

    logger.info("MAFFT workflow finished")


if __name__ == "__main__":
    main()
