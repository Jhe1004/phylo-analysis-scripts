import importlib.util
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

TRINITY_EXECUTABLE_NAME = "Trinity"
TRINITY_THREADS = "40"
TRINITY_MEMORY = "40G"
TRINITY_OUTPUT_SUFFIX = "_trinity"

RUN_TRINITY_ASSEMBLY = True
RUN_EXTRACT_LONGEST_ISOFORM = True

INPUT_FORWARD_SUFFIXES = [
    "_1.fq.gz",
    "_1.fastq.gz",
    "_1.fq",
    "_1.fastq",
    "_1.fa",
    "_1.fasta",
    "_1.fa.gz",
    "_1.fasta.gz",
]

LONGEST_INPUT_FILE_PATTERN = "*/Trinity.fasta"
LONGEST_OUTPUT_EXTENSION = ".longest.fas"
GENE_ID_REGEX = r"^(.*)_i\d+$"

DRY_RUN = False
# ============================================================


def load_helper_module(module_name, module_path):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load helper module: {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def validate_paths():
    helper_dir = os.path.join(SCRIPT_DIR, "dependencies", "scripts")
    batch_helper = os.path.join(helper_dir, "batch_trinity_helper.py")
    longest_helper = os.path.join(helper_dir, "extract_longest_isoform_helper.py")

    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not os.path.isfile(batch_helper):
        raise FileNotFoundError(f"Missing helper script: {batch_helper}")
    if not os.path.isfile(longest_helper):
        raise FileNotFoundError(f"Missing helper script: {longest_helper}")

    return batch_helper, longest_helper


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    assembly_output_dir = os.path.join(output_dir, "assemblies")
    longest_output_dir = os.path.join(output_dir, "longest_isoforms")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(assembly_output_dir, exist_ok=True)
    os.makedirs(longest_output_dir, exist_ok=True)
    logger = utils.setup_logger("TrinityPipeline", os.path.join(output_dir, "trinity_pipeline.log"))

    logger.info("Starting Trinity pipeline")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Default conda environment: {CONDA_ENV_NAME}")

    try:
        batch_helper_path, longest_helper_path = validate_paths()
        trinity_executable = resolve_executable(TRINITY_EXECUTABLE_NAME, logger)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        sys.exit(1)

    batch_helper = load_helper_module("batch_trinity_helper", batch_helper_path)
    longest_helper = load_helper_module("extract_longest_isoform_helper", longest_helper_path)

    if RUN_TRINITY_ASSEMBLY:
        logger.info("Step 1/2: Trinity batch assembly")
        batch_helper.run_batch_trinity(
            input_directory=input_dir,
            output_directory=assembly_output_dir,
            trinity_executable=trinity_executable,
            threads=TRINITY_THREADS,
            memory=TRINITY_MEMORY,
            output_suffix=TRINITY_OUTPUT_SUFFIX,
            forward_suffixes=INPUT_FORWARD_SUFFIXES,
            dry_run=DRY_RUN,
            logger=logger,
        )
    else:
        logger.info("Step 1/2 skipped: Trinity batch assembly")

    if RUN_EXTRACT_LONGEST_ISOFORM:
        logger.info("Step 2/2: Extract longest isoform")
        longest_helper.run_extract_longest_isoform(
            input_directory=assembly_output_dir,
            output_directory=longest_output_dir,
            input_file_pattern=LONGEST_INPUT_FILE_PATTERN,
            output_extension=LONGEST_OUTPUT_EXTENSION,
            gene_id_regex=GENE_ID_REGEX,
            logger=logger,
        )
    else:
        logger.info("Step 2/2 skipped: Extract longest isoform")

    logger.info("Trinity pipeline finished")


if __name__ == "__main__":
    main()
