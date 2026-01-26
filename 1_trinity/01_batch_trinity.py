import os
import sys
import argparse
import sys

# Dynamic path adjustment to find 'lib' if executed from subdirectory
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

try:
    from lib import utils
except ImportError:
    print("Error: Could not import 'lib'. Please run this script from the project root or ensure 'lib' is in PYTHONPATH.")
    sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(description="Batch Trinity Assembly")
    
    parser.add_argument("--input", "-i", default=".", help="Input directory (default: current dir)")
    parser.add_argument("--trinity-cmd", default="Trinity", help="Trinity command (default: Trinity)")
    parser.add_argument("--threads", "-t", default="40", help="Number of CPU threads (default: 40)")
    parser.add_argument("--memory", "-m", default="40G", help="Max memory (default: 40G)")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")

    return parser.parse_args()

def main():
    args = parse_args()
    input_dir = os.path.abspath(args.input)
    logger = utils.setup_logger("Trinity", os.path.join(input_dir, "trinity_batch.log"))
    
    logger.info("Starting Batch Trinity Workflow")
    logger.info(f"Input Directory: {input_dir}")
    
    # Try multiple extensions
    extensions = ['.fq.gz', '.fastq.gz', '.fq', '.fastq', '.fa', '.fasta', '.fa.gz', '.fasta.gz']
    
    samples = {}
    
    # Custom finding logic for multiple extensions since utils.find_paired_files is strictly single ext
    # This logic is specific to Trinity script original logic
    # But to use utils, we can iterate over potential suffix pairs.
    
    # Simplification: Assume user provides standardized input or we scan for common ones.
    # Let's keep the original "scan multiple" logic but adapt it.
    
    # Attempting to find pairs using utils for each ext
    found_samples = {}
    for ext in extensions:
        try:
            pairs = utils.find_paired_files(input_dir, f"_1{ext}", f"_2{ext}")
            for p in pairs:
                if p not in found_samples:
                    found_samples[p] = ext
        except Exception:
            pass # Directory might not exist or other error, handled later
            
    if not found_samples:
        logger.warning("No paired samples found.")
        sys.exit(0)

    logger.info(f"Found {len(found_samples)} samples: {list(found_samples.keys())}")
    
    for sample, ext in found_samples.items():
        left_file = os.path.join(input_dir, f"{sample}_1{ext}")
        right_file = os.path.join(input_dir, f"{sample}_2{ext}")
        output_dir = os.path.join(input_dir, f"{sample}_trinity")
        
        seq_type = "fq"
        if ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']:
            seq_type = "fa"
            
        cmd = [
            args.trinity_cmd,
            "--seqType", seq_type,
            "--left", left_file,
            "--right", right_file,
            "--CPU", args.threads,
            "--max_memory", args.memory,
            "--output", output_dir,
            "--full_cleanup"
        ]
        
        try:
           utils.run_command(cmd, logger, dry_run=args.dry_run)
           logger.info(f"Finished assembly for {sample}")
        except Exception as e:
            logger.error(f"Failed assembly for {sample}: {e}")

if __name__ == "__main__":
    main()