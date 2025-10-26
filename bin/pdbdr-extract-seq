#!/usr/bin/env python3
import os
import sys
import requests
from multiprocessing import Pool
import pdbdr as dr

def process_single_file(args):
    (filename, input_dir, pdb_suffix, cache_dir) = args
    full_path = os.path.join(input_dir, filename)
    print(f"filename: {full_path}")

    log_output = []

    # Robust PDB ID extraction: remove extension first, then split on underscore
    base = os.path.splitext(os.path.basename(filename))[0]  # e.g., "1abc_A" or "1abc"
    pdb_id = base.split("_")[0]  # "1abc"
    cached_pdb = os.path.join(cache_dir, f"{pdb_id.upper()}.pdb")

    # Figure out chain from the input file (use full path!)
    try:
        _, _, _, _, chain_of_interest = dr.pdb_to_tensor(full_path)
    except Exception as e:
        log_output.append(f"Error: Exception reading {filename} to detect chain: {e}\n")
        return (filename, "", "".join(log_output), "error", "")

    # Download the full PDB if not cached
    if not os.path.exists(cached_pdb):
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        print(url)
        try:
            response = requests.get(url)
            if response.status_code != 200:
                log_output.append(f"Error: Failed to download PDB {pdb_id} (status {response.status_code}).\n")
                return (filename, "", "".join(log_output), "error", chain_of_interest)
            with open(cached_pdb, "w") as f:
                f.write(response.text)
        except Exception as e:
            log_output.append(f"Error: Exception while downloading {pdb_id}: {e}\n")
            return (filename, "", "".join(log_output), "error", chain_of_interest)

    # Extract sequence
    try:
        _, extracted_sequence, masked_sequence, _, _ = dr.pdb_to_tensor(cached_pdb, chain_of_interest)
        print(extracted_sequence)
    except Exception as e:
        log_output.append(f"Error: Exception processing {pdb_id} using pdb_to_tensor: {e}\n")
        return (filename, "", "".join(log_output), "error", chain_of_interest)

    is_incomplete = '_' in masked_sequence
    has_noncanonical = any(c not in dr.VALID_AA for c in extracted_sequence)

    if has_noncanonical:
        category = "corrupt"
        log_output.append(f"Warning: {pdb_id} contains noncanonical characters.\n")
    elif is_incomplete:
        category = "incomplete"
    else:
        category = "complete"

    return (filename, extracted_sequence, "".join(log_output), category, chain_of_interest)

def oldprocess_single_file(args):
    (filename, input_dir, pdb_suffix, cache_dir) = args
    print(f"filename: {filename}")
    _,_,_,_,chain_of_interest = dr.pdb_to_tensor(filename)
    log_output = []
    pdb_id = filename.split("_")[0]
    cached_pdb = os.path.join(cache_dir, f"{pdb_id.upper()}.pdb")
    
    # Download the full PDB file if not cached
    if not os.path.exists(cached_pdb):
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        print(url)
        try:
            response = requests.get(url)
            if response.status_code != 200:
                log_output.append(f"Error: Failed to download PDB {pdb_id} (status {response.status_code}).\n")
                return (filename, "", "".join(log_output), "error")
            with open(cached_pdb, "w") as f:
                f.write(response.text)
        except Exception as e:
            log_output.append(f"Error: Exception while downloading {pdb_id}: {e}\n")
            return (filename, "", "".join(log_output), "error")
    
    # Extract sequence using pdb_to_tensor
    try:
        _, extracted_sequence, masked_sequence, _, _ = dr.pdb_to_tensor(cached_pdb,chain_of_interest)
        print(extracted_sequence)
    except Exception as e:
        log_output.append(f"Error: Exception processing {pdb_id} using pdb_to_tensor: {e}\n")
        return (filename, "", "".join(log_output), "error")
    is_incomplete = '_' in masked_sequence
    has_noncanonical = any(c not in dr.VALID_AA for c in extracted_sequence)
    
    # Determine output category
    if has_noncanonical:
        category = "corrupt"
        log_output.append(f"Warning: {pdb_id} contains noncanonical characters.\n")
    elif is_incomplete:
        category = "incomplete"
    else:
        category = "complete"

    return (filename, extracted_sequence, "".join(log_output), category, chain_of_interest)

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_sequence.py <input_directory> [pdb_suffix]")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)
    
    pdb_suffix = sys.argv[2] if len(sys.argv) == 3 else ".pdb"
    cache_dir = "pdb_cache"
    os.makedirs(cache_dir, exist_ok=True)
    
    complete_tsv = "complete.tsv"
    incomplete_tsv = "incomplete.tsv"
    corrupt_tsv = "corrupt.tsv"
    log_file = "extraction.log"

    #files_to_process = [f for f in os.listdir(input_dir) if f.endswith(".pdb") and pdb_suffix in f]
    # Expect pdb_suffix to include the extension when provided (e.g., ".pdb", "_A.pdb")
    pdb_suffix_norm = pdb_suffix.lower()
    files_to_process = [
        f for f in os.listdir(input_dir)
        if f.lower().endswith(pdb_suffix_norm)
    ]

    with open(log_file, "w") as log_f, \
         open(complete_tsv, "w") as comp_f, \
         open(incomplete_tsv, "w") as incomp_f, \
         open(corrupt_tsv, "w") as corrupt_f:
        
        header = "PDB_File\tChain_of_Interest\tMost_True_Sequence\n"
        comp_f.write(header)
        incomp_f.write(header)
        corrupt_f.write(header)
        
        with Pool() as pool:
            results = pool.map(
                process_single_file,
                [(f, input_dir, pdb_suffix, cache_dir) for f in files_to_process]
            )
        
        for (filename, full_sequence, log_text, category, chain_of_interest) in results:
            if category == "error":
                log_f.write(log_text)
                continue
            elif category == "complete":
                comp_f.write(f"{filename}\t{chain_of_interest}\t{full_sequence}\n")
            elif category == "incomplete":
                incomp_f.write(f"{filename}\t{chain_of_interest}\t{full_sequence}\n")
            elif category == "corrupt":
                corrupt_f.write(f"{filename}\t{chain_of_interest}\t{full_sequence}\n")

            log_f.write(log_text)

    print(f"Processing complete.\n"
          f"Complete sequences: {complete_tsv}\n"
          f"Incomplete sequences: {incomplete_tsv}\n"
          f"Corrupt sequences: {corrupt_tsv}\n"
          f"Log file: {log_file}")

if __name__ == "__main__":
    main()

