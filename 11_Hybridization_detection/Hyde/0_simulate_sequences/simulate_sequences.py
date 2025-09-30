#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
import argparse
import json
import multiprocessing
import pyvolve
from io import StringIO
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_raxml_info(log_file_path):
    """
    Parses GTR+GAMMA model parameters from a RAxML info or log file.
    Raises ValueError if any parameter is not found.
    """
    with open(log_file_path, 'r', encoding='utf-8') as f:
        log_content = f.read()

    params = {}
    alpha_match = re.search(r"alpha:\s+([\d\.]+)", log_content)
    if alpha_match:
        params['alpha'] = float(alpha_match.group(1))
    else:
        raise ValueError("Could not find 'alpha' parameter in the log file.")

    rate_names = ["A <-> C", "A <-> G", "A <-> T", "C <-> G", "C <-> T", "G <-> T"]
    rates = {}
    for name in rate_names:
        match = re.search(r"rate\s+" + name + r":\s+([\d\.]+)", log_content)
        if match:
            rates[name] = float(match.group(1))
        else:
            raise ValueError("Could not find rate parameter 'rate {}' in the log file.".format(name))
    params['rates'] = [
        rates["A <-> C"], rates["A <-> G"], rates["A <-> T"],
        rates["C <-> G"], rates["C <-> T"], rates["G <-> T"]
    ]

    freq_names = ["pi(A)", "pi(C)", "pi(G)", "pi(T)"]
    freqs = {}
    for name in freq_names:
        match = re.search(r"freq\s+" + re.escape(name) + r":\s+([\d\.]+)", log_content)
        if match:
            base = name[3]
            freqs[base] = float(match.group(1))
        else:
            raise ValueError("Could not find frequency parameter 'freq {}' in the log file.".format(name))

    if len(freqs) == 4:
        params['freqs'] = [freqs['A'], freqs['C'], freqs['G'], freqs['T']]
    else:
        raise ValueError("Failed to parse all four base frequencies.")

    return params

def generate_param_template(template_filename):
    """
    Generates a parameter template file when parsing fails.
    """
    template_data = {
        "alpha": 0.5,
        "state_freqs": [0.25, 0.25, 0.25, 0.25],
        "rates": {
            "AC": 1.0, "AG": 1.0, "AT": 1.0,
            "CG": 1.0, "CT": 1.0, "GT": 1.0
        }
    }
    with open(template_filename, 'w', encoding='utf-8') as f:
        json.dump(template_data, f, indent=4)
    print("\n[!] Parameter template file generated: '{}'".format(template_filename))
    print("[!] Please edit the parameters in this file manually and then rerun the script using the -p option.")

def parse_json_params(json_file_path):
    """
    Parses parameters from a user-provided JSON file.
    """
    with open(json_file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    if 'alpha' not in data or 'state_freqs' not in data or 'rates' not in data:
        raise ValueError("Invalid custom parameter file format. Missing 'alpha', 'state_freqs', or 'rates' key.")
    
    r = data['rates']
    rate_list = [r['AC'], r['AG'], r['AT'], r['CG'], r['CT'], r['GT']]

    return {
        'alpha': data['alpha'],
        'freqs': data['state_freqs'],
        'rates': rate_list
    }

def concatenate_fasta_files(file_list, output_file):
    """
    Concatenates multiple FASTA files into a supermatrix.
    """
    print("\n[*] Concatenating all simulated FASTA files...")
    master_sequences = {}
    all_taxa = set()
    file_lengths = {}

    for f in file_list:
        try:
            records = list(SeqIO.parse(f, "fasta"))
            if not records:
                print("    [!] Warning: File '{}' is empty and will be skipped.".format(f))
                continue
            file_lengths[f] = len(records[0].seq)
            for record in records:
                all_taxa.add(record.id)
        except FileNotFoundError:
            print("    [!] Warning: Temporary file '{}' not found. It might have failed during simulation.".format(f))


    for taxon in sorted(list(all_taxa)):
        concatenated_seq = ""
        for f in file_list:
            if f not in file_lengths: continue
            
            record_found = False
            try:
                for record in SeqIO.parse(f, "fasta"):
                    if record.id == taxon:
                        concatenated_seq += str(record.seq)
                        record_found = True
                        break
            except FileNotFoundError:
                record_found = False # File does not exist, so record cannot be found

            if not record_found:
                concatenated_seq += "?" * file_lengths[f]
        master_sequences[taxon] = Seq(concatenated_seq)

    final_records = [SeqRecord(seq, id=name, description="") for name, seq in master_sequences.items()]
    # **[MODIFIED]** Changed "fasta" to "fasta-2line" to ensure each sequence is on a single line.
    SeqIO.write(final_records, output_file, "fasta-2line")
    print("    [+] Concatenation complete. Supermatrix saved to: {}".format(output_file))

def run_simulation_worker(task_info):
    """
    A single simulation task for parallel processing.
    """
    replicate_num, model_config, tree_str, partition_config, temp_base_name = task_info
    
    try:
        model = pyvolve.Model("nucleotide", params=model_config)
        tree = pyvolve.read_tree(tree=tree_str)
        partition = pyvolve.Partition(models=model, **partition_config)
        evolver = pyvolve.Evolver(partitions=partition, tree=tree)
        
        output_filename = "{}_{}.fasta".format(temp_base_name, replicate_num)
        evolver(seqfile=output_filename, seqformat='fasta')
        
        # Be quiet in subprocesses, just return the filename
        return output_filename
    except Exception as e:
        # Catch any possible error and report it
        return "ERROR: Replicate {} failed with error: {}".format(replicate_num, e)


def main():
    """
    Main execution function.
    """
    parser = argparse.ArgumentParser(
        description="A tool to automatically simulate and concatenate sequences based on RAxML results or a custom parameter file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-i", "--info", help="Path to the info or log file generated by RAxML (automatic mode).")
    input_group.add_argument("-p", "--params", help="Path to the custom parameter file in JSON format (manual mode).")
    
    parser.add_argument("-t", "--tree", required=True, help="Path to the file containing a tree in Newick format.")
    parser.add_argument("-r", "--replicates", type=int, default=10, help="Number of simulation replicates (default: 10).")
    parser.add_argument("-c", "--chunk_size", type=int, default=100000, help="Sequence length for each simulation chunk (default: 100000).")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Number of parallel processes to run (default: 1).\nSet higher on multi-core machines for acceleration.")
    
    args = parser.parse_args()

    tree_file = args.tree
    if not os.path.exists(tree_file):
        print("\nError: Tree file '{}' not found.".format(tree_file))
        sys.exit(1)

    print("--- Sequence Simulation Script Initializing ---")
    
    model_params = None
    if args.info:
        log_file = args.info
        if not os.path.exists(log_file):
            print("\nError: Log file '{}' not found.".format(log_file))
            sys.exit(1)
        print("[*] Mode: Automatic (using log file: {})".format(log_file))
        print("[*] Using tree file: {}".format(tree_file))
        
        try:
            print("\n[*] Parsing GTR+GAMMA model parameters from log file...")
            model_params = parse_raxml_info(log_file)
            print("[+] Parameters parsed successfully!")
        except ValueError as e:
            print("\n[!] Parameter parsing failed: {}".format(e))
            generate_param_template('params_template.json')
            sys.exit(1)

    elif args.params:
        param_file = args.params
        if not os.path.exists(param_file):
            print("\nError: Parameter file '{}' not found.".format(param_file))
            sys.exit(1)
        print("[*] Mode: Manual (using parameter file: {})".format(param_file))
        print("[*] Using tree file: {}".format(tree_file))

        try:
            print("\n[*] Reading parameters from JSON file...")
            model_params = parse_json_params(param_file)
            print("[+] Parameters read successfully!")
        except (ValueError, json.JSONDecodeError) as e:
            print("\n[!] Failed to read custom parameter file: {}".format(e))
            sys.exit(1)

    r = model_params['rates']
    rate_matrix = [
        ["-", r[0], r[1], r[2]],
        [r[0], "-", r[3], r[4]],
        [r[1], r[3], "-", r[5]],
        [r[2], r[4], r[5], "-"]
    ]
    gtr_params = {"state_freqs": model_params['freqs'], "rate_matrix": rate_matrix}

    print("\n[*] Preparing for chunk-based simulation...")
    total_len = args.replicates * args.chunk_size
    print("[*] Will perform {} simulation(s), each {} bp long.".format(args.replicates, args.chunk_size))
    print("[*] Total final sequence length will be approximately: {} bp.".format(total_len))

    with open(tree_file, 'r', encoding='utf-8') as f:
        tree_str = f.read().strip()
    if not tree_str:
        print("Error: Tree file '{}' is empty or could not be read.".format(tree_file))
        sys.exit(1)
    
    partition_config = {
        'size': args.chunk_size,
        'num_rate_cats': 4,
        'gamma_shape': model_params['alpha']
    }
    temp_base_name = "__temp_sim"
    tasks = [
        (i + 1, gtr_params, tree_str, partition_config, temp_base_name)
        for i in range(args.replicates)
    ]

    print("[*] Starting simulation with {} parallel process(es)...".format(args.jobs))
    
    # Use context manager for Pool
    with multiprocessing.Pool(processes=args.jobs) as pool:
        results = pool.map(run_simulation_worker, tasks)

    print("[+] All simulation chunks completed!")

    temp_files = [res for res in results if not res.startswith("ERROR:")]
    errors = [res for res in results if res.startswith("ERROR:")]
    if errors:
        print("\n[!] The following errors occurred during simulation:")
        for err in errors:
            print("    - {}".format(err))
    
    if not temp_files:
        print("\n[!] No simulation files were generated successfully. Exiting.")
        sys.exit(1)

    base_name = os.path.splitext(os.path.basename(tree_file))[0]
    output_filename = 'simulated_supermatrix_from_{}.fasta'.format(base_name)
    
    concatenate_fasta_files(temp_files, output_filename)

    print("\n[*] Cleaning up temporary files...")
    for f in temp_files:
        try:
            os.remove(f)
        except FileNotFoundError:
            pass # Skip if file was not created due to an error
    print("    [+] Cleanup complete.")
    
    print("\n--- Task Finished ---")

if __name__ == '__main__':
    main()