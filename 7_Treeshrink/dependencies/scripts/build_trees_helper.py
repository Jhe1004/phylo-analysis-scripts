import os
import subprocess
import multiprocessing


def extract_number_from_filename(filename):
    import re

    patterns_to_try = [r"ortho(\d+)", r"(\d+)\.fas", r"(\d+)\.fasta"]
    for pattern in patterns_to_try:
        match = re.search(pattern, os.path.basename(filename))
        if match:
            return int(match.group(1))
    return float("inf")


def list_alignment_files(input_dir, pattern):
    import glob

    search_path = os.path.join(input_dir, pattern)
    files = glob.glob(search_path)
    return sorted(files, key=extract_number_from_filename)


def run_single_raxml(
    alignment_path,
    output_dir,
    raxml_executable,
    threads_per_process,
    bootstraps,
    model,
    seed,
    additional_parameters,
):
    base_name = os.path.splitext(os.path.basename(alignment_path))[0]
    command = [
        raxml_executable,
        "-T",
        str(threads_per_process),
        "-n",
        base_name,
        "-s",
        alignment_path,
        "-m",
        model,
        "-p",
        str(seed),
        "-f",
        "a",
        "-N",
        str(bootstraps),
        "-x",
        str(seed),
        *additional_parameters,
    ]

    process = subprocess.run(
        command,
        cwd=output_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=True,
    )

    tree_path = os.path.join(output_dir, f"RAxML_bestTree.{base_name}")
    return {
        "alignment": alignment_path,
        "tree": tree_path,
        "stdout": process.stdout.strip(),
        "stderr": process.stderr.strip(),
    }


def run_single_raxml_task(task):
    return run_single_raxml(**task)


def build_gene_trees(
    input_dir,
    alignment_pattern,
    output_dir,
    raxml_executable,
    process_count,
    threads_per_process,
    bootstraps,
    model,
    seed,
    additional_parameters,
    logger,
    dry_run,
):
    alignment_files = list_alignment_files(input_dir, alignment_pattern)
    if not alignment_files:
        logger.warning(f"No alignment files found with pattern: {alignment_pattern}")
        return []

    logger.info(f"Found {len(alignment_files)} alignment files for tree building")

    if dry_run:
        for alignment_path in alignment_files:
            logger.info(f"[DRY RUN] Would build tree for: {alignment_path}")
        return [os.path.join(output_dir, f"RAxML_bestTree.{os.path.splitext(os.path.basename(path))[0]}") for path in alignment_files]

    tasks = [
        {
            "alignment_path": path,
            "output_dir": output_dir,
            "raxml_executable": raxml_executable,
            "threads_per_process": threads_per_process,
            "bootstraps": bootstraps,
            "model": model,
            "seed": seed,
            "additional_parameters": additional_parameters,
        }
        for path in alignment_files
    ]

    if process_count == 1:
        results = [run_single_raxml_task(task) for task in tasks]
    else:
        with multiprocessing.Pool(processes=min(process_count, len(tasks))) as pool:
            results = pool.map(run_single_raxml_task, tasks)

    tree_paths = []
    for result in results:
        logger.info(f"Built tree: {result['tree']}")
        if result["stdout"]:
            logger.info(f"STDOUT:\n{result['stdout']}")
        if result["stderr"]:
            logger.info(f"STDERR:\n{result['stderr']}")
        tree_paths.append(result["tree"])
    return tree_paths
