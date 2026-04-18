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
INPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/0_trimmomatic/input"
OUTPUT_DIRECTORY = "/home/hejian2/My_work/disk6/phylo-analysis-scripts/0_trimmomatic/output"
CONDA_ENV_NAME = "trinity_env"

TRIMMOMATIC_JAR = os.path.join(SCRIPT_DIR, "dependencies", "trimmomatic-0.40.jar")
ADAPTERS_DIRECTORY = os.path.join(SCRIPT_DIR, "dependencies", "adapters")

FORWARD_SUFFIX = "1.fq.gz"
REVERSE_SUFFIX = "2.fq.gz"

THREADS = "30"
PHRED = "-phred33"

# If you need adapter trimming, replace this list with something like:
# TRIM_STEPS = [
#     f"ILLUMINACLIP:{os.path.join(ADAPTERS_DIRECTORY, 'TruSeq3-PE.fa')}:2:30:10",
#     "LEADING:3",
#     "TRAILING:3",
#     "SLIDINGWINDOW:4:15",
#     "HEADCROP:8",
#     "MINLEN:36",
# ]
TRIM_STEPS = [
    "LEADING:3",
    "TRAILING:3",
    "SLIDINGWINDOW:4:15",
    "HEADCROP:8",
    "MINLEN:36",
]

DRY_RUN = False
# ============================================================


def validate_config(input_dir, logger):
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    if not os.path.isfile(TRIMMOMATIC_JAR):
        raise FileNotFoundError(f"Trimmomatic jar not found: {TRIMMOMATIC_JAR}")

    if not os.path.isdir(ADAPTERS_DIRECTORY):
        logger.warning(f"Adapters directory not found: {ADAPTERS_DIRECTORY}")


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


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    os.makedirs(output_dir, exist_ok=True)
    logger = utils.setup_logger("Trimmomatic", os.path.join(output_dir, "trimmomatic.log"))

    logger.info("Starting Trimmomatic PE workflow")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Default conda environment: {CONDA_ENV_NAME}")
    logger.info(f"Trimmomatic jar: {TRIMMOMATIC_JAR}")

    try:
        validate_config(input_dir, logger)
        java_executable = resolve_executable("java", logger)
        samples = utils.find_paired_files(input_dir, FORWARD_SUFFIX, REVERSE_SUFFIX)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        sys.exit(1)

    if not samples:
        logger.warning("No paired files found.")
        sys.exit(0)

    logger.info(f"Found {len(samples)} samples: {samples}")

    for sample in samples:
        input_f = os.path.join(input_dir, f"{sample}{FORWARD_SUFFIX}")
        input_r = os.path.join(input_dir, f"{sample}{REVERSE_SUFFIX}")

        out_pair_f = os.path.join(output_dir, f"{sample}_1_paired.fq")
        out_unpair_f = os.path.join(output_dir, f"{sample}_1_unpaired.fq")
        out_pair_r = os.path.join(output_dir, f"{sample}_2_paired.fq")
        out_unpair_r = os.path.join(output_dir, f"{sample}_2_unpaired.fq")

        final_out_f = os.path.join(output_dir, f"{sample}_1.fq")
        final_out_r = os.path.join(output_dir, f"{sample}_2.fq")

        cmd = [
            java_executable,
            "-jar",
            TRIMMOMATIC_JAR,
            "PE",
            "-threads",
            THREADS,
            PHRED,
            input_f,
            input_r,
            out_pair_f,
            out_unpair_f,
            out_pair_r,
            out_unpair_r,
            *TRIM_STEPS,
        ]

        try:
            utils.run_command(cmd, logger, cwd=SCRIPT_DIR, dry_run=DRY_RUN)

            if not DRY_RUN:
                if os.path.exists(out_unpair_f):
                    os.remove(out_unpair_f)
                if os.path.exists(out_unpair_r):
                    os.remove(out_unpair_r)
                if os.path.exists(out_pair_f):
                    os.replace(out_pair_f, final_out_f)
                if os.path.exists(out_pair_r):
                    os.replace(out_pair_r, final_out_r)
                logger.info(f"Finished processing {sample}")
        except Exception as exc:
            logger.error(f"Failed to process {sample}: {exc}")


if __name__ == "__main__":
    main()
