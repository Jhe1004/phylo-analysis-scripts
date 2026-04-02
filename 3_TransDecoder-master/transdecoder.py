import multiprocessing
import os
import shutil
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
sys.path.append(PROJECT_ROOT)

try:
    from lib import utils
except ImportError:
    print("Error: Could not import 'lib'. Please ensure the project structure is intact.")
    sys.exit(1)


# ======================= CONFIGURATION =======================
INPUT_DIRECTORY = "input"
OUTPUT_DIRECTORY = "output"
CONDA_ENV_NAME = "trinity_env"

INPUT_EXTENSIONS = [".fasta", ".fa"]
WORK_SUBDIRECTORY = "work"

LONGORFS_EXECUTABLE_NAME = "TransDecoder.LongOrfs"
PREDICT_EXECUTABLE_NAME = "TransDecoder.Predict"

THREADS = max(1, multiprocessing.cpu_count())
MIN_PROTEIN_LENGTH = "100"
NO_REFINE_STARTS = True

LONGORFS_ADDITIONAL_PARAMETERS = []
PREDICT_ADDITIONAL_PARAMETERS = []

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


def find_input_files(input_dir):
    matched_files = []
    for file_name in sorted(os.listdir(input_dir)):
        if any(file_name.endswith(ext) for ext in INPUT_EXTENSIONS):
            matched_files.append(file_name)
    return matched_files


def build_sample_workdir(work_root, input_file_name):
    sample_name = input_file_name
    for ext in INPUT_EXTENSIONS:
        if sample_name.endswith(ext):
            sample_name = sample_name[: -len(ext)]
            break
    return os.path.join(work_root, sample_name)


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    work_root = os.path.join(output_dir, WORK_SUBDIRECTORY)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(work_root, exist_ok=True)
    logger = utils.setup_logger("TransDecoder", os.path.join(output_dir, "transdecoder.log"))

    logger.info("Starting TransDecoder workflow")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Default conda environment: {CONDA_ENV_NAME}")

    try:
        validate_config(input_dir)
        longorfs_executable = resolve_executable(LONGORFS_EXECUTABLE_NAME, logger)
        predict_executable = resolve_executable(PREDICT_EXECUTABLE_NAME, logger)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        sys.exit(1)

    input_files = find_input_files(input_dir)
    if not input_files:
        logger.warning(f"No input files found with extensions: {INPUT_EXTENSIONS}")
        sys.exit(0)

    logger.info(f"Found {len(input_files)} files to process: {input_files}")

    for input_file_name in input_files:
        source_path = os.path.join(input_dir, input_file_name)
        sample_workdir = build_sample_workdir(work_root, input_file_name)
        os.makedirs(sample_workdir, exist_ok=True)

        working_fasta_path = os.path.join(sample_workdir, input_file_name)
        if not os.path.exists(working_fasta_path):
            shutil.copy2(source_path, working_fasta_path)

        longorfs_cmd = [
            longorfs_executable,
            "-t",
            input_file_name,
            "-m",
            MIN_PROTEIN_LENGTH,
            "-T",
            str(THREADS),
            *LONGORFS_ADDITIONAL_PARAMETERS,
        ]

        predict_cmd = [
            predict_executable,
            "-t",
            input_file_name,
            *PREDICT_ADDITIONAL_PARAMETERS,
        ]

        if NO_REFINE_STARTS:
            predict_cmd.insert(1, "--no_refine_starts")

        try:
            logger.info(f"Running TransDecoder.LongOrfs for {input_file_name}")
            utils.run_command(longorfs_cmd, logger, cwd=sample_workdir, dry_run=DRY_RUN)
            logger.info(f"Running TransDecoder.Predict for {input_file_name}")
            utils.run_command(predict_cmd, logger, cwd=sample_workdir, dry_run=DRY_RUN)
            logger.info(f"Finished processing {input_file_name}")
        except Exception as exc:
            logger.error(f"Failed to process {input_file_name}: {exc}")


if __name__ == "__main__":
    main()
