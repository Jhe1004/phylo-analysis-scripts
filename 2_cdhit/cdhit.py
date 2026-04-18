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
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/2_cdhit/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/2_cdhit/output"
CONDA_ENV_NAME = "trinity_env"

CDHIT_EXECUTABLE_NAME = "cd-hit-est"
INPUT_EXTENSION = ".fas"
OUTPUT_SUFFIX = ".cdhit"

SEQUENCE_IDENTITY = "0.95"
THREADS = max(1, multiprocessing.cpu_count())

ADDITIONAL_PARAMETERS = [
    "-d",
    "0",
]

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
    return sorted(
        file_name
        for file_name in os.listdir(input_dir)
        if file_name.endswith(INPUT_EXTENSION)
    )


def build_output_prefix(output_dir, input_file_name):
    base_name = os.path.splitext(input_file_name)[0]
    return os.path.join(output_dir, f"{base_name}{OUTPUT_SUFFIX}")


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    os.makedirs(output_dir, exist_ok=True)
    logger = utils.setup_logger("CD-HIT", os.path.join(output_dir, "cdhit.log"))

    logger.info("Starting CD-HIT workflow")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Default conda environment: {CONDA_ENV_NAME}")

    try:
        validate_config(input_dir)
        cdhit_executable = resolve_executable(CDHIT_EXECUTABLE_NAME, logger)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        sys.exit(1)

    input_files = find_input_files(input_dir)
    if not input_files:
        logger.warning(f"No input files found with extension {INPUT_EXTENSION}")
        sys.exit(0)

    logger.info(f"Found {len(input_files)} files to process: {input_files}")

    for input_file_name in input_files:
        input_path = os.path.join(input_dir, input_file_name)
        output_prefix = build_output_prefix(output_dir, input_file_name)

        cmd = [
            cdhit_executable,
            "-i",
            input_path,
            "-o",
            output_prefix,
            "-c",
            SEQUENCE_IDENTITY,
            "-T",
            str(THREADS),
            *ADDITIONAL_PARAMETERS,
        ]

        try:
            utils.run_command(cmd, logger, cwd=SCRIPT_DIR, dry_run=DRY_RUN)
            logger.info(f"Finished processing {input_file_name}")
        except Exception as exc:
            logger.error(f"Failed to process {input_file_name}: {exc}")


if __name__ == "__main__":
    main()
