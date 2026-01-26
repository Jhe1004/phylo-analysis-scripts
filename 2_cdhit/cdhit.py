import os
import sys
import argparse
import multiprocessing
import sys

# Dynamic path adjustment to find 'lib' if executed from subdirectory
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

try:
    from lib import utils
except ImportError:
    print("Error: Could not import 'lib'. Please run this script from the project root or ensure 'lib' is in PYTHONPATH.")
    sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(description="Batch CD-HIT Clustering")
    
    parser.add_argument("--input", "-i", default=".", help="Input directory containing .fas files (default: current dir)")
    parser.add_argument("--cdhit-cmd", default="cd-hit-est", help="CD-HIT command (default: cd-hit-est)")
    parser.add_argument("--threshold", "-c", default="0.95", help="Similarity threshold (default: 0.95)")
    parser.add_argument("--threads", "-t", type=int, default=multiprocessing.cpu_count(), help="Number of parallel processes (default: all cores)")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")

    return parser.parse_args()

def process_file(file_info):
    """
    Worker function for multiprocessing.
    file_info is a tuple: (input_path, output_path, cmd_base, dry_run)
    """
    input_path, output_path, cmd_base, dry_run = file_info
    
    # Re-setup logger inside process or just handle print/exceptions?
    # For multiprocessing, it's often easier to just rely on return codes or simple prints, 
    # but let's try to use independent loggers or just let exceptions expand.
    # To keep it simple and safe with multiprocessing:
    
    cmd = cmd_base + ["-i", input_path, "-o", output_path]
    print(f"Processing: {os.path.basename(input_path)}")
    
    if dry_run:
        print(f"[DRY RUN] {' '.join(cmd)}")
        return True, "Dry Run"
    
    try:
        # We don't pass the main logger to avoid pickling issues, use a local simple one or just run_command with a dummy
        # Creating a temporary logger or just direct subprocess
        process = utils.subprocess.run(
            cmd,
            stdout=utils.subprocess.PIPE,
            stderr=utils.subprocess.PIPE,
            universal_newlines=True,
            check=True
        )
        return True, f"Success"
    except utils.subprocess.CalledProcessError as e:
        return False, f"Failed: {e.stderr}"
    except Exception as e:
        return False, f"Error: {e}"

def main():
    args = parse_args()
    input_dir = os.path.abspath(args.input)
    logger = utils.setup_logger("CD-HIT", os.path.join(input_dir, "cdhit_batch.log"))
    
    logger.info("Starting Batch CD-HIT Workflow")
    logger.info(f"Input Directory: {input_dir}")
    
    try:
        files = utils.find_files(input_dir, ".fas")
    except FileNotFoundError:
        logger.error(f"Directory not found: {input_dir}")
        sys.exit(1)
        
    if not files:
        logger.warning(f"No .fas files found in {input_dir}")
        sys.exit(0)
    
    logger.info(f"Found {len(files)} files to process.")
    
    # Prepare tasks
    tasks = []
    base_cmd = [args.cdhit_cmd, "-c", args.threshold]
    
    for f in files:
        input_path = os.path.join(input_dir, f)
        output_path = os.path.join(input_dir, f"{f}ta") # Original logic was filename + "ta"
        tasks.append((input_path, output_path, base_cmd, args.dry_run))
    
    if args.threads == 1:
        for task in tasks:
            success, msg = process_file(task)
            if success:
                logger.info(msg)
            else:
                logger.error(msg)
    else:
        logger.info(f"Using {args.threads} processes.")
        with multiprocessing.Pool(args.threads) as pool:
            results = pool.map(process_file, tasks)
            
            for i, (success, msg) in enumerate(results):
                fname = files[i]
                if success:
                    logger.info(f"{fname}: {msg}")
                else:
                    logger.error(f"{fname}: {msg}")

if __name__ == "__main__":
    main()