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

INPUT_EXTENSION = ".pep"
OUTPUT_PREFIX = "proteinortho_result"

PROTEINORTHO_EXECUTABLE_NAMES = ["proteinortho", "proteinortho6.pl"]
DIAMOND_EXECUTABLE_NAME = "diamond"
BLASTP_EXECUTABLE_NAME = "blastp"
MAKEBLASTDB_EXECUTABLE_NAME = "makeblastdb"

CPUS = os.cpu_count()
THREADS_PER_PROCESS = "1"
PREFER_DIAMOND = True

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


def resolve_executable(executable_names, logger):
    if isinstance(executable_names, str):
        executable_names = [executable_names]

    for executable_name in executable_names:
        local_candidate = os.path.join(SCRIPT_DIR, "dependencies", "bin", executable_name)
        if os.path.isfile(local_candidate) and os.access(local_candidate, os.X_OK):
            logger.info(f"Using local dependency executable: {local_candidate}")
            return local_candidate

    for executable_name in executable_names:
        conda_candidate = find_executable_in_conda_env(executable_name, CONDA_ENV_NAME)
        if conda_candidate:
            logger.info(f"Using executable from conda env '{CONDA_ENV_NAME}': {conda_candidate}")
            return conda_candidate

    for executable_name in executable_names:
        path_candidate = shutil.which(executable_name)
        if path_candidate:
            logger.info(f"Using executable from PATH: {path_candidate}")
            return path_candidate

    raise FileNotFoundError(
        f"Required executable not found: {executable_names}. "
        f"Searched in dependencies/bin, conda env '{CONDA_ENV_NAME}', and PATH."
    )


def list_input_files(input_dir):
    return sorted(file_name for file_name in os.listdir(input_dir) if file_name.endswith(INPUT_EXTENSION))


def resolve_input_dir(logger):
    preferred_input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    if os.path.isdir(preferred_input_dir):
        return preferred_input_dir

    legacy_input_files = list_input_files(SCRIPT_DIR)
    if legacy_input_files:
        logger.warning(
            f"Input directory '{preferred_input_dir}' not found. "
            f"Falling back to the script directory for backward compatibility: {SCRIPT_DIR}"
        )
        return SCRIPT_DIR

    return preferred_input_dir


def clean_pep_files(input_dir, output_dir, pep_files, logger):
    cleaned_files = []

    for pep_file in pep_files:
        input_path = os.path.join(input_dir, pep_file)
        output_path = os.path.join(output_dir, pep_file)

        with open(input_path, "r", encoding="utf-8") as f_in, open(output_path, "w", encoding="utf-8") as f_out:
            for line in f_in:
                if line.startswith(">"):
                    f_out.write(line)
                else:
                    f_out.write(line.replace("*", ""))

        cleaned_files.append(pep_file)
        logger.info(f"Cleaned peptide file: {output_path}")

    return cleaned_files


def main():
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    os.makedirs(output_dir, exist_ok=True)
    logger = utils.setup_logger("ProteinOrtho", os.path.join(output_dir, "proteinortho.log"))
    input_dir = resolve_input_dir(logger)

    logger.info("Starting ProteinOrtho workflow")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Default conda environment: {CONDA_ENV_NAME}")

    try:
        validate_config(input_dir)
        proteinortho_executable = resolve_executable(PROTEINORTHO_EXECUTABLE_NAMES, logger)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        sys.exit(1)

    input_files = list_input_files(input_dir)
    if not input_files:
        logger.warning(f"No input files found with extension {INPUT_EXTENSION}")
        sys.exit(0)

    cleaned_files = clean_pep_files(input_dir, output_dir, input_files, logger)

    search_mode = None
    binpath_dir = None

    if PREFER_DIAMOND:
        try:
            diamond_executable = resolve_executable(DIAMOND_EXECUTABLE_NAME, logger)
            search_mode = "diamond"
            binpath_dir = os.path.dirname(diamond_executable)
        except FileNotFoundError:
            logger.warning("diamond not found, fallback to blastp+")

    if search_mode is None:
        blastp_executable = resolve_executable(BLASTP_EXECUTABLE_NAME, logger)
        makeblastdb_executable = resolve_executable(MAKEBLASTDB_EXECUTABLE_NAME, logger)
        search_mode = "blastp+"
        binpath_dir = os.path.dirname(blastp_executable)
        logger.info(f"Using makeblastdb executable: {makeblastdb_executable}")

    cmd = [
        proteinortho_executable,
        f"-project={OUTPUT_PREFIX}",
        f"-p={search_mode}",
        f"-binpath={binpath_dir}",
        f"-cpus={CPUS}",
        f"-threads_per_process={THREADS_PER_PROCESS}",
        *ADDITIONAL_PARAMETERS,
        *cleaned_files,
    ]

    try:
        utils.run_command(cmd, logger, cwd=output_dir, dry_run=DRY_RUN)
        logger.info("ProteinOrtho workflow finished")
    except Exception as exc:
        logger.error(f"ProteinOrtho workflow failed: {exc}")


if __name__ == "__main__":
    main()
