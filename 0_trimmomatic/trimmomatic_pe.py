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
    parser = argparse.ArgumentParser(description="Batch Trimmomatic PE Wrapper")
    
    parser.add_argument("--input", "-i", default=".", help="Input directory (default: current dir)")
    parser.add_argument("--jar", required=True, help="Path to trimmomatic.jar")
    parser.add_argument("--forward-suffix", "-f", default="1.fq.gz", help="Forward read suffix (default: 1.fq.gz)")
    parser.add_argument("--reverse-suffix", "-r", default="2.fq.gz", help="Reverse read suffix (default: 2.fq.gz)")
    parser.add_argument("--threads", "-t", default="30", help="Number of threads (default: 30)")
    parser.add_argument("--phred", default="-phred33", choices=["-phred33", "-phred64"], help="Quality encoding (default: -phred33)")
    parser.add_argument("--trim-params", default="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36", help="Trimmomatic parameters")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")

    return parser.parse_args()

def main():
    args = parse_args()
    input_dir = os.path.abspath(args.input)
    logger = utils.setup_logger("Trimmomatic", os.path.join(input_dir, "trimmomatic.log"))

    logger.info("Starting Trimmomatic PE Workflow")
    logger.info(f"Input Directory: {input_dir}")
    
    try:
        samples = utils.find_paired_files(input_dir, args.forward_suffix, args.reverse_suffix)
    except FileNotFoundError as e:
        logger.error(e)
        sys.exit(1)

    if not samples:
        logger.warning("No paired files found.")
        sys.exit(0)
    
    logger.info(f"Found {len(samples)} samples: {samples}")

    for sample in samples:
        input_f = os.path.join(input_dir, f"{sample}{args.forward_suffix}")
        input_r = os.path.join(input_dir, f"{sample}{args.reverse_suffix}")
        
        out_pair_f = os.path.join(input_dir, f"{sample}_1_paired.fq")
        out_unpair_f = os.path.join(input_dir, f"{sample}_1_unpaired.fq")
        out_pair_r = os.path.join(input_dir, f"{sample}_2_paired.fq")
        out_unpair_r = os.path.join(input_dir, f"{sample}_2_unpaired.fq")
        
        final_out_f = os.path.join(input_dir, f"{sample}_1.fq")
        final_out_r = os.path.join(input_dir, f"{sample}_2.fq")

        cmd = [
            "java", "-jar", args.jar, "PE",
            "-threads", args.threads,
            args.phred,
            input_f, input_r,
            out_pair_f, out_unpair_f,
            out_pair_r, out_unpair_r
        ] + args.trim_params.split()
        
        try:
            utils.run_command(cmd, logger, cwd=input_dir, dry_run=args.dry_run)
            
            if not args.dry_run:
                # Cleanup and rename
                if os.path.exists(out_unpair_f): os.remove(out_unpair_f)
                if os.path.exists(out_unpair_r): os.remove(out_unpair_r)
                if os.path.exists(out_pair_f): os.rename(out_pair_f, final_out_f)
                if os.path.exists(out_pair_r): os.rename(out_pair_r, final_out_r)
                logger.info(f"Finished processing {sample}")
                
        except Exception as e:
            logger.error(f"Failed to process {sample}: {e}")

if __name__ == "__main__":
    main()