import os
import sys
import logging
import subprocess
import shutil

def setup_logger(name, log_file=None, level=logging.INFO):
    """
    Sets up a logger with the given name and log file.
    If log_file is provided, logs will be written to it.
    Logs are always printed to the console.
    """
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Console Handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File Handler
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        
    return logger

def find_paired_files(directory, f_ext, r_ext):
    """
    Finds paired files in a directory based on forward and reverse extensions.
    Returns a sorted list of base names.
    """
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")

    files = os.listdir(directory)
    base_names = set()
    
    for f in files:
        if f.endswith(f_ext):
            base_name = f[:-len(f_ext)]
            expected_r = base_name + r_ext
            if expected_r in files:
                base_names.add(base_name)
    
    return sorted(list(base_names))

def find_files(directory, ext):
    """
    Finds all files with a specific extension in a directory.
    Returns a sorted list of full filenames.
    """
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
        
    return sorted([f for f in os.listdir(directory) if f.endswith(ext)])

def run_command(cmd, logger, cwd=None, dry_run=False, check=True):
    """
    Runs a shell command.
    """
    cmd_str = " ".join(cmd) if isinstance(cmd, list) else cmd
    logger.info(f"Command: {cmd_str}")
    
    if dry_run:
        logger.info("[DRY RUN] Command skipped.")
        return 0
        
    try:
        process = subprocess.run(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=check
        )
        if process.stdout:
            logger.info(f"STDOUT:\n{process.stdout.strip()}")
        if process.stderr:
            logger.info(f"STDERR:\n{process.stderr.strip()}")
        return process.returncode
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with return code {e.returncode}")
        if e.stdout:
            logger.error(f"STDOUT:\n{e.stdout.strip()}")
        if e.stderr:
            logger.error(f"STDERR:\n{e.stderr.strip()}")
        if check:
            raise
        return e.returncode
    except Exception as e:
        logger.error(f"Execution error: {e}")
        raise

def check_dependency(tool_name):
    """
    Checks if a tool is available in the system PATH.
    """
    return shutil.which(tool_name) is not None
