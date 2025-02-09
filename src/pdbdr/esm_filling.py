"""
Module: esm_filling
Author: [Your Name]
Date: [YYYY-MM-DD]
Description:
    This module provides functions for processing protein structures
    using ESM models. It converts PDB files to tensor representations,
    fills missing structural information via either a local model or the
    Forge API, and performs sequence/structure corrections.

    The design emphasizes robust error handling, rate limiting, and clear
    documentation. Its style reflects over 20 years of expertise in
    structural bioinformatics.

    Use the unified function fill_structure(mode, pdb_id, output_pdb) to
    toggle between local and Forge backends.
    
Functions:
    pdb_to_tensor(file_path)
    fill_struct(model, data_in, max_retries=5, rate_limit=5)
    fill_structure(mode, pdb_id, output_pdb, max_retries=5, rate_limit=5)
"""

import numpy as np
import torch
import time
import random
import os
import sys
import warnings

from .constants import (
    FORGE_URL, FORGE_MDL,
    ATOM_TYPES, ATOM_ORDER, ATOM_TYPE_NUM,
    THREE_TO_ONE
)

# Suppress specific warnings about missing metadata.
warnings.filterwarnings(
    "ignore",
    message="Entity ID not found in metadata, using None as default"
)


def pdb_to_tensor(file_path):
    """
    Convert a PDB file to a tensor representation and extract sequences.

    This function parses a PDB file to identify resolved and missing
    residues, then builds a tensor of shape
    (num_residues, ATOM_TYPE_NUM, 3) containing atom coordinates.
    It also produces an extracted sequence (using THREE_TO_ONE mapping)
    and a masked sequence (missing residues are marked with "_").

    Args:
        file_path (str): Path to the PDB file.

    Returns:
        tuple: (
            torch.Tensor of coordinates,
            str: extracted one-letter sequence,
            str: masked sequence,
            dict: mapping from sequence index to residue ID
        )
    """
    # Initialize variables.
    chain_id = None
    missing_residues = {}
    seq_position_map = {}

    # Read the entire PDB file.
    with open(file_path, "r") as f:
        pdb_lines = f.readlines()

    # Determine the chain ID from the first ATOM record.
    for line in pdb_lines:
        if line.startswith("ATOM"):
            chain_id = line[21]
            break

    # If a chain is found, parse missing residues and expected sequence.
    if chain_id:
        from .pdb_io import parse_missing_residues, parse_seqres_records
        missing_residues, seq_position_map = parse_missing_residues(
            pdb_lines, chain_id
        )
        _ = parse_seqres_records(pdb_lines)  # Expected residues; unused here.

    # Collect resolved residues and their coordinates.
    resolved_residues = {}
    for line in pdb_lines:
        if not line.startswith("ATOM"):
            continue

        atom_name = line[12:16].strip()
        res_name = line[17:20].strip().upper()
        alt_loc = line[16:17].strip()
        try:
            res_seq = int(line[22:26].strip())
        except ValueError:
            continue
        insertion_code = line[26].strip()
        full_res_id = (res_seq, insertion_code)
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue

        # Skip alternate locations other than blank or "A".
        if alt_loc and alt_loc != "A":
            continue

        if full_res_id not in resolved_residues:
            resolved_residues[full_res_id] = {
                "res_name": res_name,
                "atoms": {}
            }
        resolved_residues[full_res_id]["atoms"][atom_name] = [x, y, z]

    # Combine missing and resolved residues.
    complete_residues = {**missing_residues, **resolved_residues}
    sorted_residue_ids = sorted(
        complete_residues.keys(), key=lambda x: (x[0], x[1])
    )

    extracted_sequence = []
    masked_sequence = []
    # Preallocate a 3D array for atom coordinates.
    structured_residues = np.full(
        (len(sorted_residue_ids), ATOM_TYPE_NUM, 3),
        np.nan, dtype=float
    )

    # Build sequences and fill coordinate tensor.
    for i, res_id in enumerate(sorted_residue_ids):
        if res_id in resolved_residues:
            res_name = resolved_residues[res_id]["res_name"]
            extracted_sequence.append(THREE_TO_ONE.get(res_name, "X"))
        else:
            # For missing residues, use the identity stored in missing_residues.
            res_name = missing_residues[res_id]
            extracted_sequence.append(res_name)

        # Mask missing residues with an underscore.
        if res_id in missing_residues:
            masked_sequence.append("_")
        else:
            masked_sequence.append(THREE_TO_ONE.get(res_name, "X"))

        # Fill in available atom coordinates.
        if res_id in resolved_residues:
            residue_atoms = resolved_residues[res_id]["atoms"]
            for aname, acoords in residue_atoms.items():
                if aname in ATOM_ORDER:
                    structured_residues[i, ATOM_ORDER[aname]] = acoords

        # Map sequence position to residue ID.
        seq_position_map[i + 1] = res_id

    return (
        torch.tensor(structured_residues, dtype=torch.float32),
        "".join(extracted_sequence),
        "".join(masked_sequence),
        seq_position_map,
    )


def fill_struct(model, data_in, max_retries=5, rate_limit=5):
    """
    Fill missing structure information using a provided ESM model.

    Given a PDB file path, this function converts the structure to a tensor,
    prepares a prompt with masked regions, and iteratively calls the model
    to generate a complete structure. It implements simple rate limiting.

    Args:
        model: The ESM model instance (local or remote).
        data_in (str): Path to the input PDB file.
        max_retries (int): Maximum number of retries in case of rate limits.
        rate_limit (int): Maximum allowed calls per minute.

    Returns:
        tuple: (filled protein object, sequence position map, full sequence)
               or None if filling fails.
    """
    from esm.sdk.api import ESMProtein, GenerationConfig

    # Convert PDB to tensor and extract sequences.
    x, full_seq, mskd_seq, seq_pos = pdb_to_tensor(data_in)
    pdb_indices = [i for i, char in enumerate(mskd_seq) if char != "_"]

    # Prepare template protein and token encodings.
    template_prot = ESMProtein(
        sequence=full_seq, coordinates=x,
        potential_sequence_of_concern=True
    )
    template_toks = model.encode(template_prot)
    prompt = model.encode(ESMProtein(sequence=full_seq))

    # Set default structure tokens.
    prompt.structure = torch.full_like(prompt.sequence, 4096)
    prompt.structure[0] = 4098
    prompt.structure[-1] = 4097

    # Overwrite positions where structure is available.
    for i in range(len(full_seq)):
        if i in pdb_indices:
            prompt.structure[i] = template_toks.structure[i]

    retries = 0
    backoff = 1
    while retries < max_retries:
        try:
            time.sleep(60 / rate_limit)
            complete_toks = model.generate(
                prompt,
                GenerationConfig(
                    track="structure",
                    num_steps=mskd_seq.count('_'),
                    temperature=0
                )
            )
            complete_prot = model.decode(complete_toks)
            return complete_prot, seq_pos, full_seq
        except Exception as e:
            if "429" in str(e):
                retries += 1
                sleep_time = backoff + random.uniform(0, 0.5)
                print(
                    f"Rate limit exceeded. Retrying... Attempt {retries} "
                    f"after {sleep_time:.2f}s."
                )
                time.sleep(sleep_time)
                backoff *= 2
            else:
                print(f"Unexpected error: {e}")
                break
    return None


def fill_structure(mode, pdb_id, output_pdb, max_retries=5, rate_limit=5):
    """
    Fill missing structure using either a local ESM model or the Forge API.

    This function downloads the PDB file from RCSB, instantiates the ESM model
    based on the specified mode, fills in missing structural information, applies
    post-processing (residue numbering corrections), and writes the resulting
    structure to a PDB file.

    Args:
        mode (str): "local" or "forge" to select the ESM model backend.
        pdb_id (str): The PDB identifier.
        output_pdb (str): Output file path for the filled structure.
        max_retries (int): Maximum number of retries for model generation.
        rate_limit (int): Maximum allowed calls per minute to the model API.

    Returns:
        str: The filled protein sequence.
    """
    import requests
    from .pdb_io import read_pdb, write_pdb

    # Download PDB file from RCSB.
    pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    temp_full_pdb = "temp_full.pdb"
    r = requests.get(pdb_url)
    if r.status_code != 200:
        print(f"Error fetching PDB file for {pdb_id} from RCSB")
        sys.exit(1)
    with open(temp_full_pdb, "w") as f:
        f.write(r.text)

    # Instantiate the ESM model based on the selected mode.
    if mode.lower() == "local":
        from esm.pretrained import load_local_model
        from esm.utils.constants.models import normalize_model_name
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model_name = "esm3-open"
        model_name = normalize_model_name(model_name)
        model = load_local_model(model_name, device=device)
        if device.type != "cpu":
            model = model.to(torch.bfloat16)
    elif mode.lower() == "forge":
        from esm.sdk.forge import ESM3ForgeInferenceClient
        # Read Forge token from file.
        with open("forge_token.txt", "r") as file:
            forge_tok = file.read().strip()
        if not forge_tok:
            print("Failed to read forge_token. Ensure forge_token.txt exists.")
            sys.exit(1)
        model = ESM3ForgeInferenceClient(
            model=FORGE_MDL, url=FORGE_URL, token=forge_tok
        )
    else:
        raise ValueError("Mode must be 'local' or 'forge'.")

    # Fill the structure using the model.
    filled_result = fill_struct(model, temp_full_pdb,
                                max_retries=max_retries, rate_limit=rate_limit)
    if filled_result is None:
        print("Failed to fill structure.")
        sys.exit(1)
    filled_prot, seq_pos, full_seq = filled_result

    # Write the initial filled structure.
    filled_prot.to_pdb(output_pdb)
    print(f"Filled structure ({mode}) written to {output_pdb}")

    # Post-process: Correct residue numbering based on seq_pos mapping.
    atoms = read_pdb(output_pdb)
    residue_groups = {}
    for atom in atoms:
        key = (atom["chain_id"], atom["res_seq"],
               atom["res_name"], atom["i_code"])
        residue_groups.setdefault(key, []).append(atom)
    sorted_residues = sorted(residue_groups.keys(), key=lambda x: (x[0], x[1]))
    corrected_residue_map = {}
    offset = 0
    for idx, residue in enumerate(sorted_residues):
        chain_id, res_seq, res_name, i_code = residue
        expected_idx = idx + 1  # 1-indexed
        if expected_idx in seq_pos:
            expected_res_seq, expected_i_code = seq_pos[expected_idx]
            if res_seq != expected_res_seq:
                offset = expected_res_seq - res_seq
            corrected_residue_map[residue] = (res_seq + offset, i_code)
    for atom in atoms:
        key = (atom["chain_id"], atom["res_seq"],
               atom["res_name"], atom["i_code"])
        if key in corrected_residue_map:
            atom["res_seq"], atom["i_code"] = corrected_residue_map[key]
    write_pdb(output_pdb, atoms)

    os.remove(temp_full_pdb)
    return filled_prot.sequence

