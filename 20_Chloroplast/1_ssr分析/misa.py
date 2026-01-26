#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
_______________________________________________________________________________

Program name: misa.py (A Python port of misa.pl)
Original Authors: Thomas Thiel, Sebastian Beier
Original Release date: 25/08/20 (version 2.1)
Python Port by: Gemini
Release date: 2025-08-25

_______________________________________________________________________________

DESCRIPTION: Tool for the identification and localization of
             (I)  perfect microsatellites as well as
             (II) compound microsatellites (two individual microsatellites,
                  disrupted by a certain number of bases)

SYNTAX:      python3 misa.py <FASTA file>

    <FASTAfile>      Single file in FASTA format containing the sequence(s).

    In order to specify the search criteria, an additional file containing
    the microsatellite search parameters is required named "misa.ini", which
    has the following structure:
      (a) Following a text string beginning with 'def', pairs of numbers are
          expected, whereas the first number defines the unit size and the
          second number the lower threshold of repeats for that specific unit.
      (b) Following a text string beginning with 'int' a single number defines
          the maximal number of bases between two adjacent microsatellites in
          order to specify the compound microsatellite type.
      (c) Following a text string beginning with 'GFF' a single string is expected
          either 'true' or 'false' to indicate with optional GFF(v3) output
          should be provided.
    Example:
      definition(unit_size,min_repeats):         1-10 2-6 3-5 4-5 5-5 6-5
      interruptions(max_difference_for_2_SSRs):   100
      GFF:                                        true

EXAMPLE: python3 misa.py seqs.fasta

_______________________________________________________________________________
"""

import sys
import re
from datetime import datetime
from collections import defaultdict

VERSION = "2.1"

def print_short_help():
    """Prints short help message by reading from this script's docstring."""
    with open(sys.argv[0], 'r') as f:
        docstring_lines = []
        in_docstring = False
        for line in f:
            if 'SYNTAX:' in line:
                in_docstring = True
            if in_docstring:
                match = re.match(r'^(?:SYNTAX:|EXAMPLE:|\s{4}.*)', line)
                if match:
                    docstring_lines.append(line.strip())
                else:
                    if len(docstring_lines) > 2: # Heuristic to stop after the main block
                        break
    print("\n".join(docstring_lines))


def print_full_help():
    """Prints the full help message from the docstring."""
    print(__doc__)


def parse_config(config_file="misa.ini"):
    """Parses the misa.ini configuration file."""
    try:
        with open(config_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        sys.exit(f"\nError: Specifications file '{config_file}' doesn't exist !\n")

    typrep = {}
    interrupt_max = 0
    gff_output = False

    for line in lines:
        line = line.strip()
        if re.match(r'^def\S*', line, re.IGNORECASE):
            # Original perl script was a bit strange here, this is a more robust interpretation
            # of "1-10 2-6 3-5..."
            parts = re.findall(r'(\d+)-(\d+)', line)
            for unit_size, min_repeats in parts:
                typrep[int(unit_size)] = int(min_repeats)
        elif re.match(r'^int\S*', line, re.IGNORECASE):
            match = re.search(r'(\d+)', line)
            if match:
                interrupt_max = int(match.group(1))
        elif re.match(r'^GFF\S*\s+true', line, re.IGNORECASE):
            gff_output = True
            
    # The original perl script used a regex that would produce {1=>10, 2=>6, ...}
    # from "1-10 2-6...". The following code block emulates that behavior if the
    # hyphenated format is not used.
    if not typrep:
        for line in lines:
            if re.match(r'^def\S*', line, re.IGNORECASE):
                content = re.split(r'\s+', line, 1)[1]
                numbers = re.findall(r'(\d+)', content)
                if len(numbers) % 2 == 0:
                    it = iter(numbers)
                    typrep = {int(k): int(v) for k, v in zip(it, it)}
                    break

    if not typrep:
        sys.exit("\nError: Could not parse SSR definitions from 'misa.ini'.\n"
                 "Expected format like: definition(unit_size,min_repeats): 1-10 2-6 3-5\n")

    return typrep, interrupt_max, gff_output


def fasta_parser(filename):
    """A simple FASTA parser that yields (id, sequence) tuples."""
    try:
        with open(filename, 'r') as f:
            header = None
            sequence = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if header:
                        yield header, ''.join(sequence)
                    header = line[1:]
                    sequence = []
                else:
                    sequence.append(line)
            if header:
                yield header, ''.join(sequence)
    except FileNotFoundError:
        sys.exit(f"\nError: FASTA file '{filename}' doesn't exist !\n\n")


def is_redundant(motif):
    """Check if a motif is a repetition of a smaller unit."""
    motif_len = len(motif)
    for i in range(1, motif_len // 2 + 1):
        if motif_len % i == 0:
            unit = motif[:i]
            repeats = motif_len // i
            if unit * repeats == motif:
                return True
    return False

def get_ssr_type_and_note_gff(motif_len):
    """Helper for GFF output formatting."""
    ssr_type = ""
    note_prefix = ""
    if motif_len == 1:
        ssr_type = "monomeric_repeat"
        note_prefix = "monomeric_repeat,"
    else:
        ssr_type = "microsatellite"
        if motif_len == 2:
            note_prefix = "dinucleotide_repeat_microsatellite_feature,"
        elif motif_len == 3:
            note_prefix = "trinucleotide_repeat_microsatellite_feature,"
        elif motif_len == 4:
            note_prefix = "tetranucleotide_repeat_microsatellite_feature,"
        else: # motif_len > 4
            note_prefix = "microsatellite,"
    return ssr_type, note_prefix

def get_canonical_motif(motif):
    """
    Finds the canonical representation of a motif by checking all
    circular permutations and their reverse complements.
    """
    # Generate circular permutations
    perms = {motif[i:] + motif[:i] for i in range(len(motif))}
    lex_min_fwd = min(perms)

    # Generate reverse complement and its permutations
    rev_comp = motif.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
    rev_perms = {rev_comp[i:] + rev_comp[:i] for i in range(len(rev_comp))}
    lex_min_rev = min(rev_perms)

    # The group name is the lexicographically smallest of the two, formatted as in the perl script
    if lex_min_fwd < lex_min_rev:
        return f"{lex_min_fwd}/{lex_min_rev}"
    else:
        return f"{lex_min_rev}/{lex_min_fwd}"

def main():
    """Main function to run the MISA tool."""
    # --- DECLARATION ---
    if len(sys.argv) < 2:
        print_short_help()
        sys.exit()

    if sys.argv[1].lower() in ['-h', '--h', '-help', '--help']:
        print_full_help()
        sys.exit()

    fasta_file = sys.argv[1]
    typrep, interrupt_max, gff_output = parse_config()
    
    loctime = datetime.now().strftime('%Y-%m-%d')

    # Open output file for standard .misa output if not in GFF mode
    out_misa = None
    if not gff_output:
        try:
            out_misa = open(f"{fasta_file}.misa", "w")
            out_misa.write("ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend\n")
        except IOError as e:
            sys.exit(f"Error: Could not open output file for writing: {e}")

    # --- CORE ---
    # Statistics variables
    max_repeats = 1
    min_repeats_overall = 1000
    count_motifs = defaultdict(int)
    # count_motifs_by_repeat will store counts like: {'AG': {10: 2, 12: 1}}
    count_motifs_by_repeat = defaultdict(lambda: defaultdict(int))
    count_class = defaultdict(int)
    number_sequences = 0
    size_sequences = 0
    ssr_containing_seqs_count = 0
    seqs_with_more_than_one_ssr = 0
    ssr_in_compound = 0

    # Process FASTA file
    for header, seq in fasta_parser(fasta_file):
        # Clean sequence
        original_seq = seq
        seq = re.sub(r'[\d\s>]', '', original_seq, flags=re.IGNORECASE).upper()
        
        # Clean ID
        if gff_output:
            seq_id = header.split()[0]
        else:
            seq_id = header.strip()
            seq_id = re.sub(r'\s+', '_', seq_id)

        number_sequences += 1
        size_sequences += len(seq)

        # List to store found SSRs for this sequence
        ssrs_found = []
        
        sorted_unit_sizes = sorted(typrep.keys())

        for motif_len in sorted_unit_sizes:
            min_reps = typrep[motif_len]
            if min_reps < min_repeats_overall:
                min_repeats_overall = min_reps
            
            # The regex looks for a motif of {motif_len}, then that same motif repeated {min_reps-1} or more times.
            search_regex = f"(([ACGT]{{{motif_len}}})(\\2{{{min_reps-1},}}))"
            
            for match in re.finditer(search_regex, seq, re.IGNORECASE):
                ssr_sequence = match.group(1).upper()
                motif = match.group(2).upper()

                if is_redundant(motif):
                    continue

                repeats = len(ssr_sequence) / motif_len
                end_pos = match.end()
                start_pos = end_pos - len(ssr_sequence) + 1

                ssrs_found.append({
                    "motif": motif,
                    "repeats": int(repeats),
                    "start": start_pos,
                    "end": end_pos,
                    "len": len(seq)
                })

                # Update stats
                count_motifs[motif] += 1
                count_motifs_by_repeat[motif][int(repeats)] += 1
                count_class[motif_len] += 1
                if repeats > max_repeats:
                    max_repeats = int(repeats)

        if not ssrs_found:
            continue

        ssr_containing_seqs_count += 1
        if len(ssrs_found) > 1:
            seqs_with_more_than_one_ssr += 1

        # Sort SSRs by start position
        ssrs_found.sort(key=lambda x: x['start'])
        
        # Open GFF file for this specific sequence if needed
        out_gff = None
        gff_file_opened = False
        
        # Process sorted SSRs to find simple and compound ones
        i = 0
        ssr_count_in_seq = 0
        while i < len(ssrs_found):
            current_ssr = ssrs_found[i]
            
            # Check for compound SSRs
            compound_members = [current_ssr]
            end_of_compound = current_ssr['end']
            
            j = i + 1
            while j < len(ssrs_found):
                next_ssr = ssrs_found[j]
                # interruption distance
                distance = next_ssr['start'] - end_of_compound
                
                if 1 <= distance <= interrupt_max:
                    compound_members.append(next_ssr)
                    end_of_compound = next_ssr['end']
                    j += 1
                # Overlapping case
                elif distance < 1:
                    compound_members.append(next_ssr)
                    end_of_compound = next_ssr['end']
                    j += 1
                else:
                    break
            
            ssr_count_in_seq += 1

            # Now, process what we found (either a simple or compound SSR)
            if len(compound_members) == 1: # --- Simple SSR ---
                ssr = compound_members[0]
                motif_len = len(ssr['motif'])
                
                if gff_output:
                    ssr_type, note_prefix = get_ssr_type_and_note_gff(motif_len)
                    ssr_seq_str = f"({ssr['motif']}){ssr['repeats']}"
                    note = f"{note_prefix}{ssr_seq_str}"
                    start = ssr['start']
                    end = ssr['end']
                    seq_len = ssr['len']
                    
                    if not gff_file_opened:
                        out_gff = open(f"{seq_id}.gff", "w")
                        out_gff.write("##gff-version 3\n")
                        out_gff.write(f"##sequence-region {seq_id} 1 {seq_len}\n")
                        out_gff.write(f"#!Date {loctime}\n")
                        out_gff.write(f"#!Type DNA\n")
                        out_gff.write(f"#!Source-version MISA {VERSION}\n")
                        out_gff.write(f"{seq_id}\tMISA\tregion\t1\t{seq_len}\t.\t.\t.\tID={seq_id}.1\n")
                        gff_file_opened = True
                    
                    out_gff.write(f"{seq_id}\tMISA\t{ssr_type}\t{start}\t{end}\t.\t.\t.\tNote={note};ID={seq_id}.{ssr_count_in_seq+1}\n")

                else: # Standard .misa output
                    ssr_type = f"p{motif_len}"
                    ssr_seq_str = f"({ssr['motif']}){ssr['repeats']}"
                    size = ssr['end'] - ssr['start'] + 1
                    out_misa.write(f"{seq_id}\t{ssr_count_in_seq}\t{ssr_type}\t{ssr_seq_str}\t{size}\t{ssr['start']}\t{ssr['end']}\n")

            else: # --- Compound SSR ---
                ssr_in_compound += len(compound_members)
                
                if gff_output:
                    seq_len = compound_members[0]['len']
                    if not gff_file_opened:
                        out_gff = open(f"{seq_id}.gff", "w")
                        out_gff.write("##gff-version 3\n")
                        out_gff.write(f"##sequence-region {seq_id} 1 {seq_len}\n")
                        out_gff.write(f"#!Date {loctime}\n")
                        out_gff.write(f"#!Type DNA\n")
                        out_gff.write(f"#!Source-version MISA {VERSION}\n")
                        out_gff.write(f"{seq_id}\tMISA\tregion\t1\t{seq_len}\t.\t.\t.\tID={seq_id}.1\n")
                        gff_file_opened = True

                    is_overlapping = any(compound_members[k+1]['start'] - compound_members[k]['end'] < 1 for k in range(len(compound_members)-1))
                    
                    if is_overlapping:
                        start = compound_members[0]['start']
                        end = compound_members[-1]['end']
                        ssr_type = "repeat_region"
                        note_parts = ["repeat_region"]
                        for ssr in compound_members:
                             _, note_prefix = get_ssr_type_and_note_gff(len(ssr['motif']))
                             ssr_seq_str = f"({ssr['motif']}){ssr['repeats']}"
                             note_parts.append(f"{note_prefix}{ssr_seq_str}")
                        note = ",".join(note_parts)
                        out_gff.write(f"{seq_id}\tMISA\t{ssr_type}\t{start}\t{end}\t.\t.\t.\tNote={note};ID={seq_id}.{ssr_count_in_seq+1}\n")

                    else: # Non-overlapping compound
                        for ssr in compound_members:
                            ssr_count_in_seq += 1 # Increment for each part of the compound SSR
                            motif_len = len(ssr['motif'])
                            ssr_type, note_prefix = get_ssr_type_and_note_gff(motif_len)
                            note_prefix = f"compound_repeat,{note_prefix}"
                            ssr_seq_str = f"({ssr['motif']}){ssr['repeats']}"
                            note = f"{note_prefix}{ssr_seq_str}"
                            out_gff.write(f"{seq_id}\tMISA\t{ssr_type}\t{ssr['start']}\t{ssr['end']}\t.\t.\t.\tNote={note};ID={seq_id}.{ssr_count_in_seq+1}\n")
                        ssr_count_in_seq -= 1 # Correct the counter as it was one logical block

                else: # Standard .misa output
                    ssr_seq_parts = []
                    is_overlapping_compound = False
                    for k in range(len(compound_members)):
                        ssr = compound_members[k]
                        if k == 0:
                            ssr_seq_parts.append(f"({ssr['motif']}){ssr['repeats']}")
                        else:
                            prev_ssr = compound_members[k-1]
                            distance = ssr['start'] - prev_ssr['end']
                            if distance > 1:
                                interruption = seq[prev_ssr['end']:ssr['start']-1].lower()
                                ssr_seq_parts.append(f"{interruption}({ssr['motif']}){ssr['repeats']}")
                            elif distance == 1: # adjacent
                                ssr_seq_parts.append(f"({ssr['motif']}){ssr['repeats']}")
                            else: # Overlapping
                                # This part is tricky to replicate exactly, misa.pl has a complex representation.
                                # A simple asterisk notation is more common.
                                is_overlapping_compound = True
                                ssr_seq_parts.append(f"({ssr['motif']}){ssr['repeats']}")
                    
                    ssr_type = 'c*' if is_overlapping_compound else 'c'
                    ssr_seq_str = "".join(ssr_seq_parts)
                    if is_overlapping_compound:
                        ssr_seq_str += "*" # Add trailing asterisk for compatibility
                    
                    start = compound_members[0]['start']
                    end = compound_members[-1]['end']
                    size = end - start + 1
                    out_misa.write(f"{seq_id}\t{ssr_count_in_seq}\t{ssr_type}\t{ssr_seq_str}\t{size}\t{start}\t{end}\n")


            # Move index past the members of the compound SSR
            i += len(compound_members)
            
        if out_gff:
            out_gff.close()


    # --- FINAL CLEANUP and STATISTICS ---
    if not gff_output and out_misa:
        out_misa.close()

    total_ssrs = sum(count_class.values())

    # Write statistics file
    with open(f"{fasta_file}.statistics", "w") as out_stats:
        # Specifications
        out_stats.write("Specifications\n==============\n\n")
        out_stats.write(f"Sequence source file: \"{fasta_file}\"\n\n")
        out_stats.write("Definement of microsatellites (unit size / minimum number of repeats):\n")
        spec_str = " ".join([f"({k}/{v})" for k, v in sorted(typrep.items())])
        out_stats.write(f"{spec_str}\n")
        if interrupt_max > 0:
            out_stats.write(f"\nMaximal number of bases interrupting 2 SSRs in a compound microsatellite:  {interrupt_max}\n")
        out_stats.write("\n\n\n")
        
        # Overview
        out_stats.write("RESULTS OF MICROSATELLITE SEARCH\n================================\n\n")
        out_stats.write(f"Total number of sequences examined:             {number_sequences}\n")
        out_stats.write(f"Total size of examined sequences (bp):          {size_sequences}\n")
        out_stats.write(f"Total number of identified SSRs:                {total_ssrs}\n")
        out_stats.write(f"Number of SSR containing sequences:             {ssr_containing_seqs_count}\n")
        out_stats.write(f"Number of sequences containing more than 1 SSR: {seqs_with_more_than_one_ssr}\n")
        out_stats.write(f"Number of SSRs present in compound formation:   {ssr_in_compound}\n\n\n")

        # Distribution to different repeat type classes
        out_stats.write("Distribution to different repeat type classes\n---------------------------------------------\n\n")
        out_stats.write("Unit size\tNumber of SSRs\n")
        for unit_size in sorted(count_class.keys()):
            out_stats.write(f"{unit_size}\t{count_class[unit_size]}\n")
        out_stats.write("\n")

        # Frequency of identified SSR motifs
        out_stats.write("Frequency of identified SSR motifs\n----------------------------------\n\n")
        header = "Repeats" + "".join([f"\t{i}" for i in range(min_repeats_overall, max_repeats + 1)]) + "\ttotal\n"
        out_stats.write(header)
        sorted_motifs = sorted(count_motifs.keys(), key=lambda m: (len(m), m))
        for motif in sorted_motifs:
            motif_len = len(motif)
            row = [motif]
            for j in range(min_repeats_overall, max_repeats + 1):
                if j < typrep[motif_len]:
                    row.append("-")
                else:
                    count = count_motifs_by_repeat[motif].get(j, "")
                    row.append(str(count))
            row.append(str(count_motifs[motif]))
            out_stats.write("\t".join(row) + "\n")
        out_stats.write("\n")
        
        # Frequency of classified repeat types (considering sequence complementary)
        out_stats.write("Frequency of classified repeat types (considering sequence complementary)\n-------------------------------------------------------------------------\n\n")
        
        canonical_counts = defaultdict(lambda: defaultdict(int))
        canonical_totals = defaultdict(int)
        
        for motif, counts_by_repeat in count_motifs_by_repeat.items():
            canonical_name = get_canonical_motif(motif)
            for repeat_num, count in counts_by_repeat.items():
                canonical_counts[canonical_name][repeat_num] += count
                canonical_totals[canonical_name] += count
        
        out_stats.write(header) # Same header as before
        sorted_canonicals = sorted(canonical_totals.keys(), key=lambda c: (len(c.split('/')[0]), c))
        
        for can_name in sorted_canonicals:
            motif_len = len(can_name.split('/')[0])
            row = [can_name]
            for j in range(min_repeats_overall, max_repeats + 1):
                if j < typrep.get(motif_len, float('inf')): # Use get for safety
                    row.append("-")
                else:
                    count = canonical_counts[can_name].get(j, "")
                    row.append(str(count))
            row.append(str(canonical_totals[can_name]))
            out_stats.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()