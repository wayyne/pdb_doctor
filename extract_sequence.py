#!/usr/bin/env python3
"""
Extract accurate full-length amino acid sequences from metadata‐scrubbed PDB files,
and split output into two files:
  - complete.tsv: For PDBs that already contain all residues.
  - incomplete.tsv: For PDBs where missing residues were inserted from REMARK 465.

Usage:
    python extract_sequence.py <input_directory>

Input:
    A directory containing PDB files with names in the format: <PDBID><SUFFIX>.pdb
    (Each file contains one protein chain. The chain ID is not in the filename,
    so it is determined by reading the ATOM records.)

Processing:
    For each file:
      1. Extract the PDB ID from the filename.
      2. Read the file’s ATOM records to determine the chain ID.
      3. Download the full PDB file from RCSB using the REST API.
      4. Parse the full PDB file:
           a. Extract resolved residues from ATOM records (for the chain).
           b. Parse REMARK 465 records (using fixed substrings) to get missing residues.
      5. Merge resolved and missing residues to produce the full, corrected sequence.
      6. If the PDB is missing residues (i.e. the full sequence differs from the resolved),
         write it to incomplete.tsv; otherwise, write it to complete.tsv.
      7. Produce a verbose log file (extraction.log) with details for each PDB.
         
Implementation Details:
    • Uses only standard Python libraries: requests, os, sys.
    • Caches downloaded full PDB files in a subdirectory to avoid redundant API calls.
"""

import os
import sys
import requests

# Mapping from three-letter to one-letter amino acid codes.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Additional codes can be added as needed.
}

def get_chain_id_from_file(filepath):
    """
    Read the input (metadata-scrubbed) PDB file and determine its (only) chain ID from the ATOM records.
    """
    chain_ids = set()
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    # In standard PDB format, the chain ID is at column 22 (index 21).
                    if len(line) >= 22:
                        chain = line[21]
                        chain_ids.add(chain if chain.strip() else "")
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return None
    if not chain_ids:
        return None
    elif len(chain_ids) == 1:
        return chain_ids.pop()
    else:
        # If multiple chain IDs are found (unexpected), choose one arbitrarily.
        print(f"Warning: Multiple chain IDs found in {filepath}. Using one arbitrarily.")
        return chain_ids.pop()

def download_full_pdb(pdb_id, cache_dir):
    """
    Download the full PDB file for the given pdb_id from RCSB (if not already cached).
    Returns a list of lines from the file.
    """
    pdb_filename = os.path.join(cache_dir, f"{pdb_id.upper()}.pdb")
    if os.path.exists(pdb_filename):
        with open(pdb_filename, 'r') as f:
            return f.readlines()
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        response = requests.get(url)
        if response.status_code != 200:
            print(f"Failed to download PDB {pdb_id} from RCSB (status code {response.status_code}).")
            return None
        pdb_data = response.text
        # Cache the downloaded file.
        with open(pdb_filename, 'w') as f:
            f.write(pdb_data)
        return pdb_data.splitlines()
    except Exception as e:
        print(f"Error downloading PDB {pdb_id}: {e}")
        return None

def parse_resolved_residues(pdb_lines, chain_id):
    """
    Parse the full PDB file's ATOM records (for the specified chain) and extract resolved residues.
    Returns a dictionary keyed by (resSeq, insertion_code) with values: (three_letter, one_letter).
    """
    resolved = {}
    last_key = None
    for line in pdb_lines:
        if line.startswith("ATOM"):
            if len(line) < 27:
                continue
            current_chain = line[21]
            if current_chain != chain_id:
                continue
            # Extract residue sequence number and insertion code.
            res_seq_str = line[22:26].strip()
            try:
                res_seq = int(res_seq_str)
            except ValueError:
                continue
            i_code = line[26].strip()  # Normalize insertion code (if blank, becomes empty string).
            key = (res_seq, i_code)
            if key == last_key:
                continue  # Skip duplicate entries.
            last_key = key
            res_name = line[17:20].strip().upper()
            one_letter = THREE_TO_ONE.get(res_name, "X")
            resolved[key] = (res_name, one_letter)
    return resolved

def parse_missing_residues(pdb_lines, chain_id):
    """
    Parse REMARK 465 lines from the full PDB file to extract missing residues for the specified chain.
    Uses fixed-width substrings.
    
    Expected format (by column positions, 0-indexed):
      - Columns 0-9: "REMARK 465"
      - Columns 15-18: residue name (3 letters)
      - Column 19: chain identifier
      - Columns 20-26: residue sequence number (may be signed and have spaces)
      - Column 26 (if present): insertion code (optional)
    
    Lines where the residue name (columns 15-18) is not in THREE_TO_ONE are skipped.
    Returns a dictionary keyed by (resSeq, insertion_code) with values: (three_letter, one_letter).
    """
    missing = {}
    for line in pdb_lines:
        if not line.startswith("REMARK 465"):
            continue
        if len(line) < 26:
            continue
        # Extract residue name from columns 15-18.
        res_name = line[15:18].strip().upper()
        if res_name not in THREE_TO_ONE:
            continue  # Likely a header line.
        # Extract chain ID from column 19.
        missing_chain = line[19:20].strip()
        if missing_chain != chain_id:
            continue
        # Extract residue sequence number from columns 20-26.
        res_seq_str = line[20:26].strip()
        if not res_seq_str:
            continue
        try:
            res_seq = int(res_seq_str)
        except ValueError:
            continue
        # Extract optional insertion code from column 26.
        insertion = ""
        if len(line) >= 27:
            insertion = line[26:27].strip()
        one_letter = THREE_TO_ONE.get(res_name, "X")
        key = (res_seq, insertion)
        missing[key] = (res_name, one_letter)
    return missing

def merge_residues(resolved, missing):
    """
    Merge the resolved and missing residue dictionaries.
    Resolved residues take precedence if the same key exists.
    Returns:
       - full_sequence: the final sequence (as a string) in one-letter code,
       - sorted_keys: a list of keys in order,
       - merged: the merged dictionary.
    """
    merged = resolved.copy()
    for key, value in missing.items():
        if key not in merged:
            merged[key] = value
    sorted_keys = sorted(merged.keys(), key=lambda x: (x[0], x[1]))
    full_sequence = "".join(merged[k][1] for k in sorted_keys)
    return full_sequence, sorted_keys, merged

def get_sequence_from_dict(residues_dict):
    """
    Build a sequence string (and list of sorted keys) from a residue dictionary.
    """
    sorted_keys = sorted(residues_dict.keys(), key=lambda x: (x[0], x[1]))
    sequence = "".join(residues_dict[k][1] for k in sorted_keys)
    return sequence, sorted_keys

def diff_sequences(resolved_dict, full_dict):
    """
    Identify keys (residues) that are in the full (merged) dictionary but not in the resolved dictionary.
    Returns a sorted list of such keys.
    """
    missing_keys = sorted(set(full_dict.keys()) - set(resolved_dict.keys()), key=lambda x: (x[0], x[1]))
    return missing_keys

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_sequence.py <input_directory> [pdb_suffix]")
        sys.exit(1)

    input_dir = sys.argv[1]
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)

    if len(sys.argv) == 3:
      pdb_suffix = sys.argv[2]
    else:
      pdb_suffix = ".pdb"

    # Define output file names.
    complete_tsv = "complete.tsv"
    incomplete_tsv = "incomplete.tsv"
    log_file = "extraction.log"
    cache_dir = "pdb_cache"
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    # Open output files and log file.
    with open(log_file, "w") as log_f, \
         open(complete_tsv, "w") as comp_f, \
         open(incomplete_tsv, "w") as incomp_f:

        # Write headers for TSV files.
        header = "PDB_File\tMost_True_Sequence\n"
        comp_f.write(header)
        incomp_f.write(header)

        for filename in os.listdir(input_dir):
            if not filename.endswith(".pdb"):
                continue
            if pdb_suffix not in filename:
                continue

            file_path = os.path.join(input_dir, filename)
            # Extract PDB ID from the filename (assuming the portion before the first underscore).
            pdb_id = filename.split("_")[0]

            # Determine chain ID from the metadata-scrubbed file.
            chain_id = get_chain_id_from_file(file_path)
            if chain_id is None:
                log_f.write("=" * 60 + "\n")
                log_f.write(f"File: {filename}\n")
                log_f.write("Error: No ATOM records found to determine chain ID.\n")
                log_f.write("=" * 60 + "\n\n")
                continue

            # Download full PDB file from RCSB.
            full_pdb_lines = download_full_pdb(pdb_id, cache_dir)
            if full_pdb_lines is None:
                log_f.write("=" * 60 + "\n")
                log_f.write(f"File: {filename}\n")
                log_f.write(f"Error: Failed to download full PDB for {pdb_id}.\n")
                log_f.write("=" * 60 + "\n\n")
                continue

            # Parse resolved residues from ATOM records.
            resolved = parse_resolved_residues(full_pdb_lines, chain_id)
            resolved_seq, resolved_keys = get_sequence_from_dict(resolved)

            # Parse missing residues from REMARK 465 records.
            missing = parse_missing_residues(full_pdb_lines, chain_id)

            # Merge resolved and missing residues.
            full_sequence, full_keys, merged_dict = merge_residues(resolved, missing)
            diff_keys = diff_sequences(resolved, merged_dict)

            # Write the output line to the appropriate file.
            if diff_keys:
                # Incomplete: missing residues were inserted.
                incomp_f.write(f"{filename}\t{full_sequence}\n")
            else:
                # Complete: input PDB already contained all residues.
                comp_f.write(f"{filename}\t{full_sequence}\n")

            # Write detailed log information.
            log_f.write("=" * 60 + "\n")
            log_f.write(f"Processing file: {filename}\n")
            log_f.write("-" * 60 + "\n")
            log_f.write(f"PDB ID: {pdb_id}\n")
            log_f.write(f"Chain ID: {chain_id if chain_id else '[blank]'}\n")
            log_f.write("\nResolved Residues (from ATOM records):\n")
            for k in resolved_keys:
                res_info = resolved[k]
                log_f.write(f"  Residue {k[0]}{k[1] if k[1] else ''}: {res_info[0]} -> {res_info[1]}\n")
            log_f.write("\nMissing Residues (from REMARK 465):\n")
            if missing:
                for k in sorted(missing.keys(), key=lambda x: (x[0], x[1])):
                    res_info = missing[k]
                    log_f.write(f"  Residue {k[0]}{k[1] if k[1] else ''}: {res_info[0]} -> {res_info[1]}\n")
            else:
                log_f.write("  None\n")
            log_f.write("\nReconstructed Full Sequence:\n")
            log_f.write(f"  {full_sequence}\n")
            log_f.write("\nDifferences (Missing residues inserted):\n")
            if diff_keys:
                for k in diff_keys:
                    if k in missing:
                        res_info = missing[k]
                        log_f.write(f"  Inserted missing residue at position {k[0]}{k[1] if k[1] else ''}: {res_info[0]} -> {res_info[1]}\n")
            else:
                log_f.write("  None (Resolved sequence is complete)\n")
            log_f.write("=" * 60 + "\n\n")

    print(f"Processing complete.\nComplete sequences file: {complete_tsv}\nIncomplete sequences file: {incomplete_tsv}\nLog file: {log_file}")

if __name__ == "__main__":
    main()

