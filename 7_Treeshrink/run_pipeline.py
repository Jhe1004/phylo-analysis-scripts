import glob
import importlib.util
import multiprocessing
import os
import random
import shutil
import subprocess
import sys

from Bio import SeqIO

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
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

ALIGNMENT_INPUT_SUBDIR = "alignments"
RAXML_OUTPUT_SUBDIR = "raxml_trees"
TREESHRINK_OUTPUT_SUBDIR = "treeshrink"
FINAL_ALIGNMENT_SUBDIR = "shrunk_alignments"

ALIGNMENT_PATTERN = "ortho*_cds_maffted.fas"
COMBINED_TREE_FILENAME = "combined_trees.tre"
TAXA_TO_REMOVE_FILENAME = "output.txt"

RAXML_EXECUTABLE_NAME = "raxmlHPC-PTHREADS"
RAXML_PROCESS_COUNT = max(1, multiprocessing.cpu_count() // 4)
RAXML_THREADS_PER_PROCESS = 4
RAXML_BOOTSTRAPS = 2
RAXML_MODEL = "GTRGAMMA"
RAXML_SEED = 12345
RAXML_ADDITIONAL_PARAMETERS = []

RUN_TREESHRINK_SCRIPT = os.path.join(SCRIPT_DIR, "dependencies", "scripts", "run_treeshrink.py")
TREESHRINK_EXTRA_ARGS = []
TREESHRINK_FORCE_OVERRIDE = True

FILTER_PROCESS_COUNT = max(1, multiprocessing.cpu_count() // 2)
VERIFY_RANDOM_CHECKS = 5
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


def validate_paths(input_dir, helper_path):
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not os.path.isfile(helper_path):
        raise FileNotFoundError(f"Missing helper script: {helper_path}")
    if not os.path.isfile(RUN_TREESHRINK_SCRIPT):
        raise FileNotFoundError(f"Missing TreeShrink script: {RUN_TREESHRINK_SCRIPT}")


def combine_tree_files(tree_files, combined_tree_file, logger):
    if not tree_files:
        raise FileNotFoundError("No tree files were generated for TreeShrink input.")

    with open(combined_tree_file, "w", encoding="utf-8") as outfile:
        for tree_file in tree_files:
            if not os.path.exists(tree_file):
                continue
            with open(tree_file, "r", encoding="utf-8") as infile:
                content = infile.read().strip()
                if content:
                    outfile.write(content + "\n")
    logger.info(f"Combined {len(tree_files)} gene trees into: {combined_tree_file}")


def run_treeshrink(python_executable, rscript_executable, combined_trees_file, output_dir, logger):
    command = [
        python_executable,
        RUN_TREESHRINK_SCRIPT,
        "-t",
        combined_trees_file,
        "-o",
        output_dir,
        *TREESHRINK_EXTRA_ARGS,
    ]
    if TREESHRINK_FORCE_OVERRIDE:
        command.append("--force")

    run_env = os.environ.copy()
    run_env["PATH"] = os.pathsep.join(
        [os.path.dirname(rscript_executable), os.path.dirname(python_executable), run_env.get("PATH", "")]
    )

    logger.info("Running TreeShrink")
    logger.info("Command: " + " ".join(command))

    if DRY_RUN:
        logger.info("[DRY RUN] TreeShrink skipped.")
        return

    process = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=True,
        env=run_env,
        cwd=SCRIPT_DIR,
    )
    if process.stdout:
        logger.info(f"TreeShrink STDOUT:\n{process.stdout.strip()}")
    if process.stderr:
        logger.info(f"TreeShrink STDERR:\n{process.stderr.strip()}")


def extract_number_from_filename(filename):
    import re

    patterns_to_try = [r"ortho(\d+)", r"(\d+)\.fas", r"(\d+)\.fasta", r"(\d+)\.tre"]
    for pattern in patterns_to_try:
        match = re.search(pattern, os.path.basename(filename))
        if match:
            return int(match.group(1))
    return float("inf")


def get_sorted_files(directory, pattern):
    search_path = os.path.join(directory, pattern)
    files = glob.glob(search_path)
    return sorted(files, key=extract_number_from_filename)


def filter_single_alignment(task):
    file_path, taxa_to_remove, output_dir = task
    original_records = list(SeqIO.parse(file_path, "fasta"))
    records_to_keep = [rec for rec in original_records if rec.id not in taxa_to_remove]
    output_filename = os.path.basename(file_path)
    output_path = os.path.join(output_dir, output_filename)
    SeqIO.write(records_to_keep, output_path, "fasta")
    return {
        "file_path": file_path,
        "output_path": output_path,
        "original_count": len(original_records),
        "removed_count_in_list": len(taxa_to_remove),
        "actually_removed": len(original_records) - len(records_to_keep),
        "actual_final_count": len(records_to_keep),
    }


def verify_filter_results(stats, num_to_check, logger):
    if not stats:
        logger.info("No files to verify.")
        return

    files_to_check = random.sample(stats, min(num_to_check, len(stats)))
    all_passed = True
    for stat in files_to_check:
        if stat["original_count"] - stat["actually_removed"] != stat["actual_final_count"]:
            all_passed = False
            logger.warning(f"Verification failed: {os.path.basename(stat['file_path'])}")
        else:
            logger.info(f"Verification passed: {os.path.basename(stat['file_path'])}")

    if all_passed:
        logger.info("All randomly checked filtered alignments passed verification.")
    else:
        logger.warning("Some randomly checked filtered alignments did not pass verification.")


def filter_alignments(alignment_dir, taxa_file, output_dir, logger):
    alignment_files = get_sorted_files(alignment_dir, ALIGNMENT_PATTERN)
    if not alignment_files:
        raise FileNotFoundError(f"No alignment files found with pattern: {ALIGNMENT_PATTERN}")
    if not os.path.exists(taxa_file):
        raise FileNotFoundError(f"TreeShrink taxa file not found: {taxa_file}")

    os.makedirs(output_dir, exist_ok=True)
    with open(taxa_file, "r", encoding="utf-8") as handle:
        lines_to_process = handle.readlines()

    if len(alignment_files) != len(lines_to_process):
        logger.warning(
            f"Alignment file count ({len(alignment_files)}) and taxa lines ({len(lines_to_process)}) do not match."
        )

    tasks = []
    for index, file_path in enumerate(alignment_files):
        if index >= len(lines_to_process):
            continue
        taxa_to_remove = set(lines_to_process[index].strip().split())
        tasks.append((file_path, taxa_to_remove, output_dir))

    if DRY_RUN:
        for task in tasks:
            logger.info(f"[DRY RUN] Would filter alignment: {task[0]}")
        return

    if FILTER_PROCESS_COUNT == 1:
        stats = [filter_single_alignment(task) for task in tasks]
    else:
        with multiprocessing.Pool(processes=min(FILTER_PROCESS_COUNT, len(tasks))) as pool:
            stats = pool.map(filter_single_alignment, tasks)

    for stat in stats:
        logger.info(f"Filtered: {stat['file_path']} -> {stat['output_path']}")
    verify_filter_results(stats, VERIFY_RANDOM_CHECKS, logger)


def main():
    input_dir = os.path.join(SCRIPT_DIR, INPUT_DIRECTORY)
    alignment_input_dir = os.path.join(input_dir, ALIGNMENT_INPUT_SUBDIR)
    output_dir = os.path.join(SCRIPT_DIR, OUTPUT_DIRECTORY)
    raxml_output_dir = os.path.join(output_dir, RAXML_OUTPUT_SUBDIR)
    treeshrink_output_dir = os.path.join(output_dir, TREESHRINK_OUTPUT_SUBDIR)
    final_alignment_dir = os.path.join(output_dir, FINAL_ALIGNMENT_SUBDIR)
    combined_tree_file = os.path.join(treeshrink_output_dir, COMBINED_TREE_FILENAME)
    taxa_to_remove_file = os.path.join(treeshrink_output_dir, TAXA_TO_REMOVE_FILENAME)
    helper_path = os.path.join(SCRIPT_DIR, "dependencies", "scripts", "build_trees_helper.py")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(raxml_output_dir, exist_ok=True)
    os.makedirs(treeshrink_output_dir, exist_ok=True)
    os.makedirs(final_alignment_dir, exist_ok=True)
    logger = utils.setup_logger("TreeShrinkPipeline", os.path.join(output_dir, "treeshrink_pipeline.log"))

    logger.info("Starting RAxML -> TreeShrink -> filter workflow")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Default conda environment: {CONDA_ENV_NAME}")
    logger.info(
        f"RAxML parallelism: process_count={RAXML_PROCESS_COUNT}, "
        f"threads_per_process={RAXML_THREADS_PER_PROCESS}"
    )
    logger.info(f"Filter parallelism: process_count={FILTER_PROCESS_COUNT}")

    try:
        validate_paths(alignment_input_dir, helper_path)
        raxml_executable = resolve_executable(RAXML_EXECUTABLE_NAME, logger)
        python_executable = resolve_executable("python", logger)
        rscript_executable = resolve_executable("Rscript", logger) if not DRY_RUN else None
    except FileNotFoundError as exc:
        logger.error(str(exc))
        sys.exit(1)

    helper_module = load_helper_module("build_trees_helper", helper_path)

    try:
        tree_files = helper_module.build_gene_trees(
            input_dir=alignment_input_dir,
            alignment_pattern=ALIGNMENT_PATTERN,
            output_dir=raxml_output_dir,
            raxml_executable=raxml_executable,
            process_count=RAXML_PROCESS_COUNT,
            threads_per_process=RAXML_THREADS_PER_PROCESS,
            bootstraps=RAXML_BOOTSTRAPS,
            model=RAXML_MODEL,
            seed=RAXML_SEED,
            additional_parameters=RAXML_ADDITIONAL_PARAMETERS,
            logger=logger,
            dry_run=DRY_RUN,
        )
        combine_tree_files(tree_files, combined_tree_file, logger)
        if DRY_RUN:
            logger.info("[DRY RUN] Pipeline stopped after tree-building stage.")
            sys.exit(0)
        run_treeshrink(python_executable, rscript_executable, combined_tree_file, treeshrink_output_dir, logger)
        filter_alignments(alignment_input_dir, taxa_to_remove_file, final_alignment_dir, logger)
    except Exception as exc:
        logger.error(f"Pipeline failed: {exc}")
        sys.exit(1)

    logger.info("Workflow finished")


if __name__ == "__main__":
    main()
