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

from .constants import (FORGE_URL, FORGE_MDL, ATOM_TYPES, ATOM_ORDER, 
                        ATOM_TYPE_NUM, THREE_TO_ONE, MODIFIED_RESIDUE_MAP)

# Suppress specific warnings about missing metadata.
warnings.filterwarnings(
    "ignore",
    message="Entity ID not found in metadata, using None as default"
)

"""
Module: esm_filling
Author: [Your Name]
Date: [YYYY-MM-DD]
Description:
    Provides functions for processing protein structures using ESM models.
    This includes converting PDB files to tensor representations (with chain
    filtering), filling missing structural information, and performing
    sequence/structure corrections.
"""

warnings.filterwarnings(
    "ignore",
    message="Entity ID not found in metadata, using None as default"
)

def get_letter(res_name):
    """Return one-letter code for a residue name."""
    if len(res_name) == 1:
        return res_name
    return THREE_TO_ONE.get(res_name, MODIFIED_RESIDUE_MAP.get(res_name, "X"))

# --- snip: top of file remains unchanged ---

def pdb_to_tensor(file_path, chain_id=None):
    """
    Convert a PDB file to a tensor representation and extract sequences.
    Only processes the specified chain (if provided). Parses missing and
    resolved residues, and builds a tensor of shape
    (num_residues, ATOM_TYPE_NUM, 3) with atom coordinates.
    
    Args:
        file_path (str): Path to the PDB file.
        chain_id (str): Chain identifier to filter records.
    
    Returns:
        tuple: (torch.Tensor, str, str, dict, str)
            - Coordinates tensor.
            - Extracted one-letter sequence.
            - Masked sequence (gaps marked as "_").
            - Mapping from sequence index to residue ID.
            - The chain_id actually used.
    """
    from .pdb_io import parse_missing_residues

    with open(file_path, "r") as f:
        pdb_lines = f.readlines()

    # If the PDB file contains multiple models, only consider the first model.
    if any(line.startswith("MODEL") for line in pdb_lines):
        model_lines = []
        in_first_model = False
        for line in pdb_lines:
            if line.startswith("MODEL"):
                in_first_model = True
                continue
            if line.startswith("ENDMDL") and in_first_model:
                break
            if in_first_model:
                model_lines.append(line)
        pdb_lines = model_lines

    # If no chain specified, derive from first ATOM record.
    if chain_id is None:
        for line in pdb_lines:
            if line.startswith("ATOM"):
                chain_id = line[21].strip()
                print(f"found chainID: {chain_id}")
                break
    if chain_id is None:
        raise ValueError("No chain ID found in PDB file and not provided.")

    missing_residues, seq_position_map = parse_missing_residues(
        pdb_lines, chain_id
    )

    resolved_residues = {}
    # Process ATOM records.
    for line in pdb_lines:
        if not line.startswith("ATOM"):
            continue
        if line[21].strip() != chain_id:
            continue
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip().upper()
        alt_loc = line[16:17].strip()
        try:
            res_seq = int(line[22:26].strip())
        except ValueError:
            continue
        insertion_code = line[26].strip()
        key = (res_seq, insertion_code)
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue
        if alt_loc and alt_loc != "A":
            continue
        if key not in resolved_residues:
            resolved_residues[key] = {"res_name": res_name, "atoms": {}}
        resolved_residues[key]["atoms"][atom_name] = [x, y, z]

    # Process HETATM records for common modified residues.
    for line in pdb_lines:
        if not line.startswith("HETATM"):
            continue
        if line[21].strip() != chain_id:
            continue
        res_name = line[17:20].strip().upper()
        if res_name not in MODIFIED_RESIDUE_MAP:
            continue
        alt_loc = line[16:17].strip()
        if alt_loc and alt_loc != "A":
            continue
        try:
            res_seq = int(line[22:26].strip())
        except ValueError:
            continue
        insertion_code = line[26].strip()
        key = (res_seq, insertion_code)
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue
        canon = MODIFIED_RESIDUE_MAP.get(res_name)
        if key not in resolved_residues:
            resolved_residues[key] = {"res_name": canon, "atoms": {}}
        atom_name = line[12:16].strip()
        resolved_residues[key]["atoms"].setdefault(atom_name, [x, y, z])

    complete_residues = {**missing_residues, **resolved_residues}
    sorted_keys = sorted(complete_residues.keys(), key=lambda x: (x[0], x[1]))

    extracted_sequence = []
    masked_sequence = []
    structured_residues = np.full((len(sorted_keys), ATOM_TYPE_NUM, 3),
                                  np.nan, dtype=float)

    for i, res_id in enumerate(sorted_keys):
        if res_id in resolved_residues:
            res_name = resolved_residues[res_id]["res_name"]
            letter = get_letter(res_name)
            extracted_sequence.append(letter)
        else:
            res_name = missing_residues[res_id]
            extracted_sequence.append(res_name)
        if res_id in missing_residues:
            masked_sequence.append("_")
        else:
            letter = get_letter(res_name)
            masked_sequence.append(letter)
        if res_id in resolved_residues:
            atoms = resolved_residues[res_id]["atoms"]
            for aname, acoords in atoms.items():
                if aname in ATOM_ORDER:
                    structured_residues[i, ATOM_ORDER[aname]] = acoords
        seq_position_map[i + 1] = res_id

    return (torch.tensor(structured_residues, dtype=torch.float32),
            "".join(extracted_sequence), "".join(masked_sequence),
            seq_position_map, chain_id)


def fill_struct(model, data_in, chain_id=None, max_retries=5, rate_limit=5):
    """
    Fill missing structure information using a provided ESM model.

    Args:
        model: The ESM model instance (local or remote).
        data_in (str): Path to the input PDB file.
        chain_id (str|None): Chain to use. If None, inferred from PDB.
        max_retries (int): Maximum number of retries in case of rate limits.
        rate_limit (int): Maximum allowed calls per minute.

    Returns:
        tuple: (filled protein object, sequence position map, full sequence)
               or None if filling fails.
    """
    from esm.sdk.api import ESMProtein, GenerationConfig

    # Convert PDB to tensor and extract sequences.
    x, full_seq, mskd_seq, seq_pos, chain_id = pdb_to_tensor(data_in, chain_id)
    pdb_indices = [i for i, char in enumerate(mskd_seq) if char != "_"]

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

def _filter_pdb_by_chain_text(pdb_text: str, chain_id: str) -> str:
    out_lines = []
    for line in pdb_text.splitlines(True):
        rec = line[0:6].strip()
        # keep ONLY ATOM/ANISOU/TER for the requested chain
        if rec in ("ATOM", "ANISOU", "TER"):
            if line[21:22] == chain_id or (chain_id == "" and line[21:22] == " "):
                out_lines.append(line)
        elif rec in ("MODEL", "ENDMDL"):
            continue
        else:
            out_lines.append(line)
    if any(l.startswith(("ATOM", "ANISOU")) for l in out_lines) and not any(l.startswith("TER") for l in out_lines):
        out_lines.append("TER\n")
    out_lines.append("END\n")
    return "".join(out_lines)


def fill_structure(mode, pdb_id, output_pdb, chain_id=None, max_retries=5, rate_limit=5):
    """
    Fill missing structure using either a local ESM model or the Forge API.

    Args:
        mode (str): "local" or "forge" to select the ESM model backend.
        pdb_id (str): The PDB identifier.
        output_pdb (str): Output file path for the filled structure.
        chain_id (str|None): Chain to keep from the downloaded PDB. If None,
            the first ATOM chain is used.
        max_retries (int): Maximum number of retries for model generation.
        rate_limit (int): Maximum allowed calls per minute to the model API.

    Returns:
        str: The filled protein sequence.
    """
    import requests
    from .pdb_io import read_pdb, write_pdb

    pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    temp_full_pdb = "temp_full.pdb"

    r = requests.get(pdb_url)
    if r.status_code != 200:
        print(f"Error fetching PDB file for {pdb_id} from RCSB")
        sys.exit(1)

    pdb_text = r.text

    # If chain_id provided, strip other chains before any processing.
    # If not provided, infer from first ATOM and then strip others.
    if chain_id is None:
        inferred = None
        for line in pdb_text.splitlines():
            if line.startswith("ATOM"):
                inferred = line[21].strip()
                break
        if inferred is None:
            print("No ATOM records found to infer chain ID.")
            sys.exit(1)
        chain_id = inferred

    filtered_text = _filter_pdb_by_chain_text(pdb_text, chain_id)
    with open(temp_full_pdb, "w") as f:
        f.write(filtered_text)

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

    # Fill the structure using the model (only the selected chain).
    filled_result = fill_struct(
        model, temp_full_pdb, chain_id=chain_id,
        max_retries=max_retries, rate_limit=rate_limit
    )
    if filled_result is None:
        print(f"Failed to fill structure: {pdb_id} chain {chain_id}.")
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
        chain_id_atom, res_seq, res_name, i_code = residue
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

