import os
import subprocess # Preferred module for running external commands
import sys      # For exiting script if needed

# --- User Configuration ---
# Specify the exact suffixes for your paired-end reads
forward_suffix = "1.fq.gz"  # Example: "_1.fastq", "_R1.fq.gz", "_forward.fastq.gz"
reverse_suffix = "2.fq.gz"  # Example: "_2.fastq", "_R2.fq.gz", "_reverse.fastq.gz"

# Define Trimmomatic output suffixes (can be customized if needed)
paired_f_out_suffix = "_1_paired.fq"
unpaired_f_out_suffix = "_1_unpaired.fq"
paired_r_out_suffix = "_2_paired.fq"
unpaired_r_out_suffix = "_2_unpaired.fq"

# Define final renamed output suffixes (can be customized if needed)
final_f_suffix = "_1.fq"
final_r_suffix = "_2.fq"

# Trimmomatic settings
trimmomatic_jar = "trimmomatic-0.39.jar" # Path to your Trimmomatic JAR file
threads = "30"
phred = "-phred33" # or "-phred64"
trim_options = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36"
# --- End User Configuration ---


def find_read_pairs(directory, f_suffix, r_suffix):
    """
    Finds pairs of sequencing files based on specified suffixes.

    Args:
        directory (str): The directory to search for files.
        f_suffix (str): The suffix of the forward read files.
        r_suffix (str): The suffix of the reverse read files.

    Returns:
        list: A sorted list of unique base names for found pairs.
              Returns an empty list if no pairs are found.
    """
    base_names = set()
    all_files = os.listdir(directory)

    print(f"Searching for files ending with '{f_suffix}' and '{r_suffix}' in '{directory}'...")

    # Create a set of all files for efficient lookup
    file_set = set(all_files)

    for filename in all_files:
        # Check if the file is a potential forward read file
        if filename.endswith(f_suffix):
            # Derive the base name by removing the forward suffix
            base_name = filename[:-len(f_suffix)]

            # Construct the expected reverse read filename
            expected_r_file = base_name + r_suffix

            # Check if the corresponding reverse read file exists
            if expected_r_file in file_set:
                base_names.add(base_name)
            else:
                print(f"  Warning: Found forward file '{filename}' but missing corresponding reverse file '{expected_r_file}'. Skipping.")

    if not base_names:
        print("  No complete pairs found.")
    else:
        print(f"  Found {len(base_names)} potential pairs.")

    # Return sorted list for consistent processing order
    return sorted(list(base_names))

# --- Main script execution ---

# Get the current working directory
current_dir = os.getcwd()

# Find the base names of the paired files
file_base_names = find_read_pairs(current_dir, forward_suffix, reverse_suffix)

# Exit if no pairs were found
if not file_base_names:
    print("\nNo matching paired-end files found to process. Exiting.")
    sys.exit(1) # Exit with a non-zero code indicating an issue

print(f"\nFound {len(file_base_names)} pairs to process: {file_base_names}")

# --- Run Trimmomatic ---
print("\n--- Running Trimmomatic ---")
trimmomatic_failed_pairs = [] # Keep track of pairs where Trimmomatic failed

for base in file_base_names:
    # Construct input filenames using defined suffixes
    input_f = os.path.join(current_dir, base + forward_suffix)
    input_r = os.path.join(current_dir, base + reverse_suffix)

    # Construct output filenames using defined output suffixes
    output_pf = os.path.join(current_dir, base + paired_f_out_suffix)
    output_uf = os.path.join(current_dir, base + unpaired_f_out_suffix)
    output_pr = os.path.join(current_dir, base + paired_r_out_suffix)
    output_ur = os.path.join(current_dir, base + unpaired_r_out_suffix)

    # Build the command as a list of arguments for subprocess
    cmd = [
        "java", "-jar", trimmomatic_jar, "PE",
        "-threads", threads, phred,
        input_f, input_r,
        output_pf, output_uf, output_pr, output_ur,
    ]
    # Add trimming options as separate arguments if they contain spaces internally,
    # but here Trimmomatic expects them as single string arguments like "SLIDINGWINDOW:4:15"
    cmd.extend(trim_options.split()) # Split options like "LEADING:3" etc.

    print(f"\nProcessing pair: {base}")
    print(f"  Command: {' '.join(cmd)}") # Print the command for clarity

    try:
        # Run the command. check=True raises CalledProcessError if command fails
        # capture_output=True captures stdout/stderr. text=True decodes them as text.
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=current_dir)
        # Trimmomatic often prints summary info to stderr, so print both
        print(f"  Trimmomatic Stdout for {base}:\n{result.stdout.strip()}")
        print(f"  Trimmomatic Stderr for {base}:\n{result.stderr.strip()}")
        print(f"  Successfully processed {base}")

    except subprocess.CalledProcessError as e:
        # Handle errors if Trimmomatic returns a non-zero exit code
        print(f"  ERROR: Trimmomatic failed for pair '{base}'!")
        print(f"  Return code: {e.returncode}")
        print(f"  Stdout: {e.stdout.strip()}")
        print(f"  Stderr: {e.stderr.strip()}")
        trimmomatic_failed_pairs.append(base) # Mark this pair as failed
    except FileNotFoundError:
        # Handle error if 'java' or the .jar file is not found
        print(f"  ERROR: Cannot find 'java' or '{trimmomatic_jar}'.")
        print("  Make sure Java is installed and in your PATH, and the JAR file path is correct.")
        print("  Exiting script.")
        sys.exit(1)
    except Exception as e:
        # Catch any other unexpected errors during subprocess execution
        print(f"  An unexpected error occurred running Trimmomatic for {base}: {e}")
        trimmomatic_failed_pairs.append(base)


# --- Cleanup and Rename ---
print("\n--- Cleaning up and Renaming Files ---")

for base in file_base_names:
    # Skip cleanup/rename if Trimmomatic failed for this pair
    if base in trimmomatic_failed_pairs:
        print(f"\nSkipping cleanup/rename for failed pair: {base}")
        continue

    print(f"\nProcessing cleanup/rename for: {base}")

    # Files to remove
    file_to_remove_uf = os.path.join(current_dir, base + unpaired_f_out_suffix)
    file_to_remove_ur = os.path.join(current_dir, base + unpaired_r_out_suffix)

    # Files to rename (source)
    file_to_rename_pf = os.path.join(current_dir, base + paired_f_out_suffix)
    file_to_rename_pr = os.path.join(current_dir, base + paired_r_out_suffix)

    # Files to rename (destination)
    final_name_f = os.path.join(current_dir, base + final_f_suffix)
    final_name_r = os.path.join(current_dir, base + final_r_suffix)

    # Remove unpaired files (with checks)
    for file_path in [file_to_remove_uf, file_to_remove_ur]:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"  Removed: {os.path.basename(file_path)}")
            else:
                print(f"  Warning: File not found for removal (already removed or not created?): {os.path.basename(file_path)}")
        except OSError as e:
            print(f"  Error removing file '{os.path.basename(file_path)}': {e}")

    # Rename paired files (with checks)
    rename_pairs = [
        (file_to_rename_pf, final_name_f),
        (file_to_rename_pr, final_name_r)
    ]
    for old_path, new_path in rename_pairs:
        try:
            if os.path.exists(old_path):
                os.rename(old_path, new_path)
                print(f"  Renamed: {os.path.basename(old_path)} -> {os.path.basename(new_path)}")
            else:
                 # It's possible Trimmomatic produced 0 paired reads output if input was very bad
                 print(f"  Warning: File not found for renaming (Trimmomatic might not have produced it?): {os.path.basename(old_path)}")
        except OSError as e:
            print(f"  Error renaming file '{os.path.basename(old_path)}' to '{os.path.basename(new_path)}': {e}")


print("\n--- Script Finished ---")