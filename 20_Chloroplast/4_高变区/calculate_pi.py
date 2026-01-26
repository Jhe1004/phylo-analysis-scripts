import sys
import argparse
from collections import Counter

def read_fasta(file_path):
    """
    Reads an aligned FASTA file and returns a dictionary of sequences.
    Handles multi-line FASTA format.
    """
    sequences = {}
    current_id = None
    current_seq = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        sequences[current_id] = ''.join(current_seq)
    return sequences

def calculate_pi_for_alignment(sequences, window_size=None, step_size=None):
    """
    Calculates nucleotide diversity (pi) for an entire alignment or in sliding windows.

    Args:
        sequences (dict): A dictionary of aligned sequences {id: seq}.
        window_size (int, optional): The size of the sliding window.
        step_size (int, optional): The step size for the sliding window.

    Returns:
        If no window parameters are given, returns the overall pi value (float).
        If window parameters are given, returns a list of tuples: (midpoint, pi_value).
    """
    if not sequences:
        raise ValueError("Input sequence dictionary is empty.")

    seq_list = list(sequences.values())
    num_sequences = len(seq_list)
    alignment_length = len(seq_list[0])

    if num_sequences < 2:
        print("Warning: Nucleotide diversity requires at least 2 sequences. Returning 0.", file=sys.stderr)
        return 0.0 if not window_size else []

    # Transpose the alignment for easy column-wise iteration
    transposed_alignment = list(zip(*seq_list))

    site_pi_values = []
    for i in range(alignment_length):
        column = transposed_alignment[i]
        
        # Filter out gaps ('-') and ambiguous bases ('N')
        valid_bases = [base.upper() for base in column if base.upper() in ('A', 'T', 'C', 'G')]
        n_k = len(valid_bases)

        if n_k < 2:
            site_pi_values.append(0.0)
            continue

        # Count base frequencies
        freqs = Counter(valid_bases)
        
        # Calculate sum of squared frequencies
        sum_p_squared = sum((count / n_k) ** 2 for count in freqs.values())
        
        # Calculate pi for the site with small sample size correction
        pi_site = (n_k / (n_k - 1)) * (1 - sum_p_squared)
        site_pi_values.append(pi_site)

    if not window_size:
        # Calculate overall pi
        if not site_pi_values:
            return 0.0
        overall_pi = sum(site_pi_values) / alignment_length
        return overall_pi
    else:
        # Calculate pi in sliding windows
        if not step_size:
            step_size = window_size
        
        window_results = []
        for start in range(0, alignment_length - window_size + 1, step_size):
            end = start + window_size
            window_slice = site_pi_values[start:end]
            
            # Check if the window is valid (some tools require a certain proportion of valid sites)
            if not window_slice:
                continue

            window_pi = sum(window_slice) / len(window_slice)
            midpoint = start + window_size // 2
            window_results.append((midpoint, window_pi))
        
        return window_results

def main():
    """Main function to parse arguments and run the calculation."""
    parser = argparse.ArgumentParser(
        description="Calculate nucleotide diversity (π) from a multiple sequence alignment in FASTA format.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "fasta_file",
        metavar="FASTA_FILE",
        type=str,
        help="Path to the aligned FASTA file."
    )
    parser.add_argument(
        "-w", "--window-size",
        type=int,
        help="Optional: The size of the sliding window for analysis. If not provided, calculates overall π."
    )
    parser.add_argument(
        "-s", "--step-size",
        type=int,
        help="Optional: The step size for the sliding window. Defaults to the window size if not provided."
    )
    
    args = parser.parse_args()

    try:
        sequences = read_fasta(args.fasta_file)
        
        if not sequences:
            print(f"Error: No sequences found in {args.fasta_file}", file=sys.stderr)
            sys.exit(1)

        # Check for consistent alignment length
        first_len = len(next(iter(sequences.values())))
        if not all(len(seq) == first_len for seq in sequences.values()):
            print("Error: Sequences in the alignment have inconsistent lengths.", file=sys.stderr)
            sys.exit(1)

        if args.window_size:
            # Sliding window analysis
            if args.window_size > first_len:
                print("Error: Window size cannot be larger than the alignment length.", file=sys.stderr)
                sys.exit(1)
            results = calculate_pi_for_alignment(sequences, args.window_size, args.step_size)
            print("Window_Midpoint\tPi_Value")
            for midpoint, pi_value in results:
                print(f"{midpoint}\t{pi_value:.8f}")
        else:
            # Overall analysis
            overall_pi = calculate_pi_for_alignment(sequences)
            print("Overall Nucleotide Diversity (π):")
            print(f"{overall_pi:.8f}")

    except FileNotFoundError:
        print(f"Error: The file '{args.fasta_file}' was not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()