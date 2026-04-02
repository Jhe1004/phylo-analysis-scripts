import os


def detect_paired_samples(input_directory, forward_suffixes):
    files = set(os.listdir(input_directory))
    found_samples = {}

    for forward_suffix in forward_suffixes:
        if not forward_suffix.startswith("_1"):
            continue
        reverse_suffix = "_2" + forward_suffix[2:]

        for filename in files:
            if filename.endswith(forward_suffix):
                sample_name = filename[: -len(forward_suffix)]
                reverse_name = sample_name + reverse_suffix
                if reverse_name in files and sample_name not in found_samples:
                    found_samples[sample_name] = (forward_suffix, reverse_suffix)

    return found_samples


def infer_seq_type(forward_suffix):
    fasta_suffixes = {".fa", ".fasta", ".fa.gz", ".fasta.gz"}
    for suffix in fasta_suffixes:
        if forward_suffix.endswith(suffix):
            return "fa"
    return "fq"


def run_batch_trinity(
    input_directory,
    output_directory,
    trinity_executable,
    threads,
    memory,
    output_suffix,
    forward_suffixes,
    dry_run,
    logger,
):
    os.makedirs(output_directory, exist_ok=True)
    found_samples = detect_paired_samples(input_directory, forward_suffixes)

    if not found_samples:
        logger.warning("No paired samples found.")
        return

    logger.info(f"Found {len(found_samples)} samples: {list(found_samples.keys())}")

    for sample_name, (forward_suffix, reverse_suffix) in found_samples.items():
        left_file = os.path.join(input_directory, f"{sample_name}{forward_suffix}")
        right_file = os.path.join(input_directory, f"{sample_name}{reverse_suffix}")
        output_dir = os.path.join(output_directory, f"{sample_name}{output_suffix}")
        seq_type = infer_seq_type(forward_suffix)

        cmd = [
            trinity_executable,
            "--seqType",
            seq_type,
            "--left",
            left_file,
            "--right",
            right_file,
            "--CPU",
            threads,
            "--max_memory",
            memory,
            "--output",
            output_dir,
            "--full_cleanup",
        ]

        logger.info(f"Running Trinity for sample: {sample_name}")
        logger.info("Command: " + " ".join(cmd))

        if dry_run:
            logger.info("[DRY RUN] Command skipped.")
            continue

        try:
            import subprocess

            process = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )
            if process.stdout:
                logger.info(f"STDOUT:\n{process.stdout.strip()}")
            if process.stderr:
                logger.info(f"STDERR:\n{process.stderr.strip()}")
            logger.info(f"Finished assembly for {sample_name}")
        except subprocess.CalledProcessError as exc:
            logger.error(f"Failed assembly for {sample_name}")
            if exc.stdout:
                logger.error(f"STDOUT:\n{exc.stdout.strip()}")
            if exc.stderr:
                logger.error(f"STDERR:\n{exc.stderr.strip()}")
