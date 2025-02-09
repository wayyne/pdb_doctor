"""
Module: pdb_io
Author: [Your Name]
Description:
    This module provides utilities to read, write, and process PDB files.
    It includes functions to parse missing residues, update chain IDs,
    group atoms by residues, and convert residue codes.
    The module is designed to be robust and maintainable, reflecting
    two decades of professional experience in structural biology data
    processing.
"""

import numpy as np
import torch
import warnings
from .constants import (
    FORGE_URL, FORGE_MDL,
    ATOM_TYPES, ATOM_ORDER, ATOM_TYPE_NUM,
    THREE_TO_ONE, ONE_TO_THREE
)

# Set print options for NumPy and PyTorch to display full arrays/tensors
np.set_printoptions(threshold=np.inf, linewidth=200, suppress=True)
torch.set_printoptions(threshold=10_000, linewidth=200)

# Suppress specific warnings about missing metadata
warnings.filterwarnings(
    "ignore",
    message="Entity ID not found in metadata, using None as default"
)

def parse_missing_residues(pdb_lines, chain_id):
    """
    Parse PDB file lines for missing or unresolved residues on a given chain.
    
    This function examines REMARK 465 records and gaps in ATOM records. For
    gaps, it computes the C-N distance between adjacent residues. A gap is
    flagged if the distance is unusually large (> 2.0 Å).
    
    Args:
        pdb_lines (list): Lines read from a PDB file.
        chain_id (str): Chain identifier to process.
    
    Returns:
        tuple:
            - missing_all (dict): Mapping from (res_seq, insertion) to one-letter
              residue code (or "?" if unknown).
            - seq_position_map (dict): Mapping from sequence index (starting at 1)
              to residue identifiers.
    """
    missing_remark = {}
    # Process REMARK 465 lines to capture missing residues
    for line in pdb_lines:
        if not line.startswith("REMARK 465"):
            continue
        if len(line) < 27:
            continue
        res_name_3 = line[15:18].strip().upper()
        if res_name_3 not in THREE_TO_ONE:
            continue
        if line[19:20].strip() != chain_id:
            continue
        res_seq_str = line[20:26].strip()
        if not res_seq_str:
            continue
        try:
            res_seq = int(res_seq_str)
        except ValueError:
            continue
        insertion = line[26:27].strip() if len(line) >= 27 else ""
        missing_remark[(res_seq, insertion)] = THREE_TO_ONE[res_name_3]

    resolved = {}
    # Process ATOM records to collect coordinates and residue info
    for line in pdb_lines:
        if not line.startswith("ATOM"):
            continue
        if len(line) < 22 or line[21].strip() != chain_id:
            continue
        res_seq_str = line[22:26].strip()
        try:
            res_seq = int(res_seq_str)
        except ValueError:
            continue
        insertion = line[26:27].strip() if len(line) >= 27 else ""
        res_name_3 = line[17:20].strip().upper()
        res_name = THREE_TO_ONE.get(res_name_3, "?")
        atom_name = line[12:16].strip()
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coord = (x, y, z)
        except ValueError:
            coord = None

        key = (res_seq, insertion)
        if key not in resolved:
            resolved[key] = {
                "res_seq": res_seq,
                "insertion": insertion,
                "res_name": res_name,
                "N": None,
                "C": None
            }
        if atom_name == "N":
            resolved[key]["N"] = coord
        elif atom_name == "C":
            resolved[key]["C"] = coord

    # Sort residues by sequence number and insertion code
    sorted_keys = sorted(resolved.keys(), key=lambda k: (k[0], k[1]))
    sorted_res = [resolved[k] for k in sorted_keys]

    missing_gap = {}

    def compute_distance(coord1, coord2):
        """
        Compute the Euclidean distance between two 3D coordinates.
        
        Returns None if either coordinate is missing.
        """
        if coord1 is None or coord2 is None:
            return None
        dx = coord1[0] - coord2[0]
        dy = coord1[1] - coord2[1]
        dz = coord1[2] - coord2[2]
        return (dx * dx + dy * dy + dz * dz) ** 0.5

    # Detect gaps by checking for non-consecutive residue numbers
    for i in range(len(sorted_res) - 1):
        res1 = sorted_res[i]
        res2 = sorted_res[i + 1]
        expected_next = res1["res_seq"] + 1
        if res2["res_seq"] > expected_next:
            for missing_seq in range(expected_next, res2["res_seq"]):
                key = (missing_seq, "")
                if key in missing_remark:
                    continue
                if res1["C"] is not None and res2["N"] is not None:
                    bond_len = compute_distance(res1["C"], res2["N"])
                    if bond_len is not None and bond_len < 2.0:
                        continue
                    else:
                        missing_gap[key] = "?"
                else:
                    missing_gap[key] = "?"

    missing_all = {}
    missing_all.update(missing_remark)
    missing_all.update(missing_gap)

    # Build a complete residue map (both missing and resolved)
    complete_residues = {**missing_all, **resolved}
    sorted_residue_ids = sorted(
        complete_residues.keys(), key=lambda x: (x[0], x[1])
    )

    seq_position_map = {}
    for i, res_id in enumerate(sorted_residue_ids, start=1):
        seq_position_map[i] = res_id

    return missing_all, seq_position_map


def insertion_order(i_code):
    """
    Compute the ordering index for an insertion code.
    
    An empty or None insertion code is assigned order 0.
    
    Args:
        i_code (str): Insertion code character.
    
    Returns:
        int: Numeric order for sorting.
    """
    if i_code == "" or i_code is None:
        return 0
    return ord(i_code.upper()) - ord("A") + 1


def read_pdb(filename):
    """
    Read a PDB file and extract atom records.
    
    Parses both ATOM and HETATM records. Alternate location indicators that
    are not blank or 'A' are skipped.
    
    Args:
        filename (str): Path to the PDB file.
    
    Returns:
        list: A list of dictionaries, each representing an atom.
    """
    atoms = []
    with open(filename, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            alt_loc = line[16].strip()
            if alt_loc not in ["", "A"]:
                continue
            try:
                record_name = line[0:6].strip()
                serial = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_seq = int(line[22:26].strip())
                i_code = line[26].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occ_str = line[54:60].strip()
                occupancy = float(occ_str) if occ_str else 1.0
                tf_str = line[60:66].strip()
                temp_factor = float(tf_str) if tf_str else 0.0
                element = line[76:78].strip()
                charge = line[78:80].strip()
            except Exception as e:
                print("Error parsing line:", line, e)
                continue
            atom = {
                "record_name": record_name,
                "serial": serial,
                "atom_name": atom_name,
                "alt_loc": alt_loc,
                "res_name": res_name,
                "chain_id": chain_id,
                "res_seq": res_seq,
                "i_code": i_code,
                "x": x,
                "y": y,
                "z": z,
                "occupancy": occupancy,
                "temp_factor": temp_factor,
                "element": element,
                "charge": charge,
            }
            atoms.append(atom)
    return atoms


def write_pdb(filename, atoms):
    """
    Write a list of atom dictionaries to a PDB file.
    
    The atoms are sorted by chain, residue sequence, insertion code,
    and atom name. Serial numbers are reassigned sequentially.
    
    Args:
        filename (str): Output file path.
        atoms (list): List of atom dictionaries.
    """
    def sort_key(atom):
        return (
            atom["chain_id"],
            atom["res_seq"],
            insertion_order(atom["i_code"]),
            atom["atom_name"],
        )

    atoms_sorted = sorted(atoms, key=sort_key)
    # Reassign serial numbers sequentially
    for i, atom in enumerate(atoms_sorted, start=1):
        atom["serial"] = i
    with open(filename, "w") as f:
        for atom in atoms_sorted:
            line = (
                "{:<6}{:>5} {:<4}{:1}{:<3} {:1}{:>4}{:<1}   "
                .format(
                    atom["record_name"],
                    atom["serial"],
                    atom["atom_name"].ljust(4),
                    atom["alt_loc"],
                    atom["res_name"],
                    atom["chain_id"],
                    atom["res_seq"],
                    atom["i_code"] if atom["i_code"] else "",
                )
            )
            line += "{:>8.3f}{:>8.3f}{:>8.3f}".format(
                atom["x"], atom["y"], atom["z"]
            )
            line += "{:>6.2f}{:>6.2f}          ".format(
                atom["occupancy"], atom["temp_factor"]
            )
            line += "{:<2}{:<2}".format(atom["element"], atom["charge"])
            f.write(line + "\n")
        f.write("END\n")


def update_chain_id(atoms, new_chain_id):
    """
    Update the chain identifier for all atoms in the list.
    
    Args:
        atoms (list): List of atom dictionaries.
        new_chain_id (str): The new chain ID to assign.
    
    Returns:
        list: The updated list of atoms.
    """
    for atom in atoms:
        atom["chain_id"] = new_chain_id
    return atoms


def group_by_residue(atoms):
    """
    Group atoms by residue.
    
    Atoms are grouped using a tuple key of (chain_id, res_seq, i_code).
    Only ATOM records are considered.
    
    Args:
        atoms (list): List of atom dictionaries.
    
    Returns:
        dict: Mapping from residue keys to lists of atoms.
    """
    residues = {}
    for atom in atoms:
        if atom["record_name"] != "ATOM":
            continue
        key = (atom["chain_id"], atom["res_seq"], atom["i_code"])
        residues.setdefault(key, []).append(atom)
    return residues


def get_backbone_atoms(res_atoms):
    """
    Extract backbone atoms from a residue.
    
    Only atoms with names 'N', 'CA', 'C', and 'O' are selected.
    
    Args:
        res_atoms (list): List of atom dictionaries for a residue.
    
    Returns:
        dict: Mapping from backbone atom names to atom dictionaries.
    """
    backbone = {}
    for atom in res_atoms:
        if atom["atom_name"] in ["N", "CA", "C", "O"]:
            backbone[atom["atom_name"]] = atom
    return backbone


def get_ca_coordinate(atom):
    """
    Retrieve the coordinate of a C-alpha atom.
    
    Args:
        atom (dict): Atom dictionary.
    
    Returns:
        np.ndarray: 3D coordinate as an array.
    """
    return np.array([atom["x"], atom["y"], atom["z"]])


def get_atom_coordinate(atom):
    """
    Retrieve the coordinate of an atom.
    
    Args:
        atom (dict): Atom dictionary.
    
    Returns:
        np.ndarray: 3D coordinate as an array.
    """
    return np.array([atom["x"], atom["y"], atom["z"]])


def parse_seqres_records(pdb_lines):
    """
    Parse SEQRES records to build the expected sequence per chain.
    
    Args:
        pdb_lines (list): Lines from a PDB file.
    
    Returns:
        dict: Mapping from chain IDs to a one-letter residue sequence.
    """
    seqres_data = {}
    for line in pdb_lines:
        if not line.startswith("SEQRES"):
            continue
        chain_id = line[11:12].strip()
        res_names = line[19:].split()
        if chain_id not in seqres_data:
            seqres_data[chain_id] = []
        seqres_data[chain_id].extend(res_names)

    # Convert three-letter codes to a one-letter sequence
    expected_seq = {}
    for chain, residues in seqres_data.items():
        expected_seq[chain] = "".join(
            THREE_TO_ONE.get(res, "X") for res in residues
        )
    return expected_seq

