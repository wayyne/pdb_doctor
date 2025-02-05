#!/usr/bin/env python3
"""
fold_and_transfer.py

This script implements a workflow to:
  1. Fold a given protein sequence (provided on the command line) using ESM3.
  2. Read an experimental (partial) PDB structure.
  3. Globally align the folded structure onto the partial structure using a custom
     Kabsch algorithm (using only CA atoms).
  4. Identify contiguous missing residue segments (per chain) and “stitch” them locally
     using selected anchor atoms (for example, using CA and O for anchors).
  5. Insert only those missing segments (backbone atoms only) into the partial structure.
  6. Perform an amide–bond sanity check.
  7. Write the globally fitted folded structure to "fitted_folded.pdb" (for inspection) and
     the final chimeric structure (with corrected atom numbering) to the specified output file.

This version now handles:
  - **Alt–locations:** Only atoms with alt_loc "" or "A" are used.
  - **Insertion codes:** Sorting and grouping now use a helper function so that residues
    with insertion codes are ordered properly.

Usage:
  python fold_and_transfer.py <partial_pdb> "<sequence>" <output_pdb>

Dependencies:
  - numpy
  - torch and the ESM3 modules (for folding)
"""

import numpy as np
import sys, os

########################################
# Helper function for insertion code ordering
########################################

def insertion_order(i_code):
    """Convert an insertion code to an integer order (blank = 0, A = 1, B = 2, etc.)."""
    if i_code == "" or i_code is None:
        return 0
    else:
        return ord(i_code.upper()) - ord('A') + 1

########################################
# Helper: Three‐letter to one‐letter conversion
########################################

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_partial_sequence(res_dict, chain):
    """
    Extracts a one-letter sequence from the residues in the given chain (from a residue dictionary).
    Returns a tuple (seq, keys) where:
      - seq is the concatenated one-letter codes (in order of increasing (res_seq, i_code))
      - keys is the corresponding ordered list of residue keys.
    """
    keys = [k for k in res_dict if k[0] == chain]
    keys = sorted(keys, key=lambda k: (k[1], insertion_order(k[2])))
    seq = ""
    for key in keys:
        res_atoms = res_dict[key]
        # Assume all atoms in a residue have the same res_name.
        res_name = res_atoms[0]["res_name"].upper()
        letter = three_to_one.get(res_name, 'X')
        seq += letter
    return seq, keys

def find_subsequence_indices(full_seq, subseq):
    """
    Finds the best alignment indices of subseq in full_seq using dynamic programming.
    Returns a list of indices corresponding to an optimal left-to-right match.
    """

    len_full = len(full_seq)
    len_sub = len(subseq)
    
    # DP table to store best alignment scores
    dp = np.zeros((len_sub + 1, len_full + 1))
    backtrack = np.full((len_sub + 1, len_full + 1), -1)  # Store backtracking info

    # Fill DP table
    for i in range(1, len_sub + 1):
        for j in range(1, len_full + 1):
            if subseq[i - 1] == full_seq[j - 1]:  # If characters match
                dp[i][j] = dp[i - 1][j - 1] + 1  # Increase match score
                backtrack[i][j] = j - 1  # Store index

            else:
                dp[i][j] = dp[i][j - 1]  # Carry forward best match
                backtrack[i][j] = backtrack[i][j - 1]

    # Backtrack to get indices
    indices = []
    j = len_full  # Start from the last column
    for i in range(len_sub, 0, -1):
        j = backtrack[i][j]
        if j == -1:
            return None  # No valid subsequence match
        indices.append(j)
    
    return indices[::-1]  # Reverse to get correct order

def oldfind_subsequence_indices(full_seq, subseq):
    """
    Given the full sequence (full_seq) and a subsequence (subseq), returns a list of indices
    (0-indexed) in full_seq that correspond to a left-to-right greedy match of subseq.
    Assumes that subseq is indeed a subsequence of full_seq.
    Returns None if the matching fails.
    """
    indices = []
    start = 0
    for char in subseq:
        pos = full_seq.find(char, start)
        if pos == -1:
            return None
        indices.append(pos)
        start = pos + 1
    return indices


########################################
# 1. PDB File I/O Functions
########################################
def update_chain_id(atoms, new_chain_id):
    """
    Updates the chain_id of all atoms in the list to new_chain_id.
    
    Parameters:
    - atoms (list of dicts): List of atom dictionaries (from read_pdb).
    - new_chain_id (str): The new chain ID to set for all atoms.

    Returns:
    - updated_atoms (list of dicts): The modified list with updated chain IDs.
    """
    for atom in atoms:
        atom["chain_id"] = new_chain_id  # Update the chain ID
    return atoms  # Return updated list (in-place modification)

def read_pdb(filename):
    """
    Read a PDB file and return a list of atom dictionaries.
    Only atoms with alt_loc "" or "A" are included.
    """
    atoms = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                alt_loc = line[16].strip()
                if alt_loc not in ["", "A"]:
                    continue  # skip alternate conformations other than primary
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
                    "charge": charge
                }
                atoms.append(atom)
    return atoms

def write_pdb(filename, atoms):
    """
    Write a list of atom dictionaries to a PDB file.
    Atoms are sorted by (chain, res_seq, insertion_order(i_code), atom_name).
    Atom serial numbers are reassigned sequentially.
    """
    def sort_key(atom):
        return (atom["chain_id"], atom["res_seq"], insertion_order(atom["i_code"]), atom["atom_name"])
    atoms_sorted = sorted(atoms, key=sort_key)
    for i, atom in enumerate(atoms_sorted, start=1):
        atom["serial"] = i
    with open(filename, "w") as f:
        for atom in atoms_sorted:
            line = "{:<6}{:>5} {:<4}{:1}{:<3} {:1}{:>4}{:<1}   ".format(
                atom["record_name"],
                atom["serial"],
                atom["atom_name"].ljust(4),
                atom["alt_loc"],
                atom["res_name"],
                atom["chain_id"],
                atom["res_seq"],
                atom["i_code"] if atom["i_code"] else ""
            )
            line += "{:>8.3f}{:>8.3f}{:>8.3f}".format(atom["x"], atom["y"], atom["z"])
            line += "{:>6.2f}{:>6.2f}          ".format(atom["occupancy"], atom["temp_factor"])
            line += "{:<2}{:<2}".format(atom["element"], atom["charge"])
            f.write(line + "\n")
        f.write("END\n")

########################################
# 2. Grouping and Backbone Extraction
########################################

def group_by_residue(atoms):
    """
    Group atoms by residue key: (chain_id, res_seq, i_code).
    Returns a dictionary mapping each key to a list of atoms.
    """
    residues = {}
    for atom in atoms:
        key = (atom["chain_id"], atom["res_seq"], atom["i_code"])
        residues.setdefault(key, []).append(atom)
    return residues

def get_backbone_atoms(res_atoms):
    """
    From a list of atoms in a residue, return a dictionary mapping
    backbone atom names (N, CA, C, O) to their atom dictionary.
    """
    backbone = {}
    for atom in res_atoms:
        if atom["atom_name"] in ["N", "CA", "C", "O"]:
            backbone[atom["atom_name"]] = atom
    return backbone

def get_ca_coordinate(atom):
    """Return the (x,y,z) coordinate of a CA atom as a NumPy array."""
    return np.array([atom["x"], atom["y"], atom["z"]])

def get_atom_coordinate(atom):
    """Return the (x,y,z) coordinate of an atom as a NumPy array."""
    return np.array([atom["x"], atom["y"], atom["z"]])

########################################
# 3. Global Alignment (Kabsch using CA atoms)
########################################

def compute_centroid(coords):
    """Compute the centroid of an array of coordinates (N x 3)."""
    return np.mean(coords, axis=0)

def kabsch(P, Q):
    """
    Compute the optimal rotation (R) and translation (t) to align P onto Q.
    P and Q are (N x 3) arrays.
    """
    centroid_P = compute_centroid(P)
    centroid_Q = compute_centroid(Q)
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    C = np.dot(P_centered.T, Q_centered)
    V, S, Wt = np.linalg.svd(C)
    d = np.linalg.det(np.dot(Wt.T, V.T))
    if d < 0:
        Wt[-1, :] *= -1
    R = np.dot(Wt.T, V.T)
    t = centroid_Q - np.dot(R, centroid_P)
    return R, t

def transform_atoms(atoms, R, t):
    """Apply transformation (R, t) to each atom (in place)."""
    for atom in atoms:
        coord = np.array([atom["x"], atom["y"], atom["z"]])
        new_coord = np.dot(R, coord) + t
        atom["x"], atom["y"], atom["z"] = new_coord.tolist()

def compute_alignment_rmsd(partial_atoms, folded_atoms, common_keys):
    """Compute RMSD over CA atoms for the residues in common_keys."""
    partial_res = group_by_residue(partial_atoms)
    folded_res = group_by_residue(folded_atoms)
    sq_errors = []
    for key in common_keys:
        pb = get_backbone_atoms(partial_res[key])
        fb = get_backbone_atoms(folded_res[key])
        if "CA" in pb and "CA" in fb:
            diff = get_ca_coordinate(pb["CA"]) - get_ca_coordinate(fb["CA"])
            sq_errors.append(np.dot(diff, diff))
    return np.sqrt(np.mean(sq_errors))

def align_folded_to_partial(partial_atoms, folded_atoms):
    """
    Identify residues (by key) that have CA atoms in both structures.
    Compute the global alignment transformation (R, t) from folded to partial.
    Return (R, t, common_keys). Raise an error if fewer than 3 common CA atoms.
    """
    partial_res = group_by_residue(partial_atoms)
    folded_res = group_by_residue(folded_atoms)
    common_keys = []
    partial_coords = []
    folded_coords = []
    for key in partial_res:
        if key in folded_res:
            pb = get_backbone_atoms(partial_res[key])
            fb = get_backbone_atoms(folded_res[key])
            if "CA" in pb and "CA" in fb:
                common_keys.append(key)
                partial_coords.append(get_ca_coordinate(pb["CA"]))
                folded_coords.append(get_ca_coordinate(fb["CA"]))
    if len(partial_coords) < 3:
        raise ValueError("Not enough common CA atoms for global alignment (need at least 3).")
    partial_coords = np.array(partial_coords)
    folded_coords = np.array(folded_coords)
    R, t = kabsch(folded_coords, partial_coords)
    return R, t, common_keys

########################################
# 4. Local Stitching Functions
########################################

def group_missing_segments(folded_keys, partial_keys):
    """
    For each chain, group contiguous folded residue keys that are missing in partial.
    Returns a dict mapping chain_id to a list of segments (each segment is a list of keys).
    """
    segments = {}
    chain_keys = {}
    for key in folded_keys:
        chain_keys.setdefault(key[0], []).append(key)
    for chain, keys in chain_keys.items():
        keys = sorted(keys, key=lambda k: (k[1], insertion_order(k[2])))
        missing = [k for k in keys if k not in partial_keys]
        segs = []
        if missing:
            current_seg = [missing[0]]
            for prev, curr in zip(missing, missing[1:]):
                if curr[1] - prev[1] == 1 and insertion_order(curr[2]) - insertion_order(prev[2]) in [0, 1]:
                    current_seg.append(curr)
                else:
                    segs.append(current_seg)
                    current_seg = [curr]
            segs.append(current_seg)
            segments[chain] = segs
    return segments

def oldlocal_stitch_terminal(segment_keys, folded_res_dict, partial_res_dict):
    """
    For a terminal gap (e.g. missing N-terminus), use the following:
      - Target anchors (from partial): CA of the first present residue and N of the following residue.
      - Source anchors (from folded): CA of the last missing residue and N of the first present residue.
    Compute the rigid transformation over these two pairs (using CA and N) and return the transformed segment and RMSD.
    """
    chain = segment_keys[0][0]
    # Determine target anchors: sort partial keys for this chain.
    partial_keys = {k for k in partial_res_dict if k[0] == chain}
    sorted_partial = sorted(list(partial_keys), key=lambda k: (k[1], insertion_order(k[2])))
    if len(sorted_partial) < 2:
        raise ValueError("Not enough partial residues for terminal gap fitting.")
    target_anchor1 = get_ca_coordinate(get_backbone_atoms(partial_res_dict[sorted_partial[0]])["CA"])
    target_anchor2 = get_atom_coordinate(get_backbone_atoms(partial_res_dict[sorted_partial[1]])["N"])
    Q = np.array([target_anchor1, target_anchor2])
    # Source anchors from folded.
    source_anchor1 = get_ca_coordinate(get_backbone_atoms(folded_res_dict[segment_keys[-1]])["CA"])
    # Find the first folded residue that is present (not missing) for this chain.
    folded_keys = sorted([k for k in folded_res_dict if k[0] == chain], key=lambda k: (k[1], insertion_order(k[2])))
    source_anchor2 = None
    for k in folded_keys:
        if k not in segment_keys:
            pb = get_backbone_atoms(folded_res_dict[k])
            if "N" in pb:
                source_anchor2 = get_atom_coordinate(pb["N"])
                break
    if source_anchor2 is None:
        raise ValueError("Could not find source anchor 2 for terminal gap.")
    P = np.array([source_anchor1, source_anchor2])
    R_local, t_local = kabsch(P, Q)
    seg_atoms = []
    for key in segment_keys:
        seg_atoms.extend(folded_res_dict[key])
    adjusted_atoms = []
    for atom in seg_atoms:
        new_atom = atom.copy()
        coord = np.array([atom["x"], atom["y"], atom["z"]])
        new_coord = np.dot(R_local, coord) + t_local
        new_atom["x"], new_atom["y"], new_atom["z"] = new_coord.tolist()
        adjusted_atoms.append(new_atom)
    transformed_points = np.dot(R_local, P.T).T + t_local
    local_rmsd = np.sqrt(np.mean(np.sum((transformed_points - Q)**2, axis=1)))
    return adjusted_atoms, local_rmsd

def local_stitch_internal(segment_keys, folded_res_dict, partial_res_dict):
    """
    For an internal gap (e.g. missing residues 100-110 with anchors at 99 and 111), use:
      - Left target anchor: from partial, use CA and C from residue 99.
      - Right target anchor: from partial, use CA from residue 111 and N from residue 112.
      - Left source anchor: from folded, use CA and C from residue 100.
      - Right source anchor: from folded, use CA from residue 110 and N from residue 111.
    Compute the rigid transformation over these four points and return the adjusted segment and RMSD.
    """
    left_source_key = segment_keys[0]
    right_source_key = segment_keys[-1]
    chain = segment_keys[0][0]
    left_target_key = (chain, segment_keys[0][1] - 1, "")
    right_target_key = (chain, segment_keys[-1][1] + 1, "")
    if left_target_key not in partial_res_dict or right_target_key not in partial_res_dict:
        raise ValueError("Missing anchor(s) in partial structure for internal gap.")
    left_target_bb = get_backbone_atoms(partial_res_dict[left_target_key])
    right_target_bb = get_backbone_atoms(partial_res_dict[right_target_key])
    sorted_partial = sorted([k for k in partial_res_dict if k[0] == chain], key=lambda k: (k[1], insertion_order(k[2])))
    right_anchor_key = None
    for idx, k in enumerate(sorted_partial):
        if k == right_target_key and idx + 1 < len(sorted_partial):
            right_anchor_key = sorted_partial[idx+1]
            break
    if right_anchor_key is None:
        raise ValueError("Could not find residue after right target anchor for internal gap.")
    right_anchor_bb = get_backbone_atoms(partial_res_dict[right_anchor_key])
    if "CA" not in left_target_bb or "C" not in left_target_bb:
        raise ValueError("Left target anchor missing CA or C.")
    target_left1 = get_ca_coordinate(left_target_bb["CA"])
    target_left2 = get_atom_coordinate(left_target_bb["C"])
    if "CA" not in right_target_bb or "N" not in right_anchor_bb:
        raise ValueError("Right target anchor missing CA or N.")
    target_right1 = get_ca_coordinate(right_target_bb["CA"])
    target_right2 = get_atom_coordinate(right_anchor_bb["N"])
    Q = np.array([target_left1, target_left2, target_right1, target_right2])
    left_source_bb = get_backbone_atoms(folded_res_dict[left_source_key])
    right_source_bb = get_backbone_atoms(folded_res_dict[right_source_key])
    if "CA" not in left_source_bb or "C" not in left_source_bb or "CA" not in right_source_bb or "N" not in right_source_bb:
        raise ValueError("Source missing required backbone atoms for internal gap.")
    source_left1 = get_ca_coordinate(left_source_bb["CA"])
    source_left2 = get_atom_coordinate(left_source_bb["C"])
    source_right1 = get_ca_coordinate(right_source_bb["CA"])
    source_right2 = get_atom_coordinate(right_source_bb["N"])
    P = np.array([source_left1, source_left2, source_right1, source_right2])
    R_local, t_local = kabsch(P, Q)
    seg_atoms = []
    for key in segment_keys:
        seg_atoms.extend(folded_res_dict[key])
    adjusted_atoms = []
    for atom in seg_atoms:
        new_atom = atom.copy()
        coord = np.array([atom["x"], atom["y"], atom["z"]])
        new_coord = np.dot(R_local, coord) + t_local
        new_atom["x"], new_atom["y"], new_atom["z"] = new_coord.tolist()
        adjusted_atoms.append(new_atom)
    transformed_points = np.dot(R_local, P.T).T + t_local
    local_rmsd = np.sqrt(np.mean(np.sum((transformed_points - Q)**2, axis=1)))
    return adjusted_atoms, local_rmsd

def local_stitch_terminal(segment_keys, folded_res_dict, partial_res_dict, terminal='N'):
    """
    For a terminal gap (either N- or C-terminal) use two anchor pairs:
      - For an N-terminal gap (terminal='N'):
          * Target anchors (from partial): CA of the first present residue and N of the second.
          * Source anchors (from folded): CA of the last missing residue and N from the first present (non‐missing) folded residue.
      - For a C-terminal gap (terminal='C'):
          * Target anchors (from partial): C of the last present residue and CA of the second-to-last.
          * Source anchors (from folded): N of the first missing residue and C from the last present (non‐missing) folded residue.
    Computes the local rigid transformation (via Kabsch) and returns the transformed missing segment and local RMSD.
    """
    chain = segment_keys[0][0]
    if terminal == 'N':
        # Use the first two residues in partial (lowest residue numbers)
        partial_keys = [k for k in partial_res_dict if k[0] == chain]
        sorted_partial = sorted(partial_keys, key=lambda k: (k[1], insertion_order(k[2])))
        if len(sorted_partial) < 2:
            raise ValueError("Not enough partial residues for N-terminal gap fitting.")
        target_anchor1 = get_ca_coordinate(get_backbone_atoms(partial_res_dict[sorted_partial[0]])["CA"])
        target_anchor2 = get_atom_coordinate(get_backbone_atoms(partial_res_dict[sorted_partial[1]])["N"])
        # Source anchors: use CA of the last missing residue and find the first folded residue (after the gap) with N.
        source_anchor1 = get_ca_coordinate(get_backbone_atoms(folded_res_dict[segment_keys[-1]])["CA"])
        folded_keys = sorted([k for k in folded_res_dict if k[0] == chain],
                             key=lambda k: (k[1], insertion_order(k[2])))
        source_anchor2 = None
        for k in folded_keys:
            if k not in segment_keys:
                pb = get_backbone_atoms(folded_res_dict[k])
                if "N" in pb:
                    source_anchor2 = get_atom_coordinate(pb["N"])
                    break
        if source_anchor2 is None:
            raise ValueError("Could not find source anchor 2 for N-terminal gap.")
        P = np.array([source_anchor1, source_anchor2])
        Q = np.array([target_anchor1, target_anchor2])
    elif terminal == 'C':
        # Use the last two residues in partial (highest residue numbers)
        partial_keys = [k for k in partial_res_dict if k[0] == chain]
        sorted_partial = sorted(partial_keys, key=lambda k: (k[1], insertion_order(k[2])))
        if len(sorted_partial) < 2:
            raise ValueError("Not enough partial residues for C-terminal gap fitting.")
        target_anchor1 = get_atom_coordinate(get_backbone_atoms(partial_res_dict[sorted_partial[-1]])["C"])
        target_anchor2 = get_ca_coordinate(get_backbone_atoms(partial_res_dict[sorted_partial[-2]])["CA"])
        # Source anchors: use N from the first missing residue and find the last folded residue (before the gap) with C.
        source_anchor1 = get_atom_coordinate(get_backbone_atoms(folded_res_dict[segment_keys[0]])["N"])
        folded_keys = sorted([k for k in folded_res_dict if k[0] == chain],
                             key=lambda k: (k[1], insertion_order(k[2])))
        source_anchor2 = None
        for k in reversed(folded_keys):
            if k not in segment_keys:
                pb = get_backbone_atoms(folded_res_dict[k])
                if "C" in pb:
                    source_anchor2 = get_atom_coordinate(pb["C"])
                    break
        if source_anchor2 is None:
            raise ValueError("Could not find source anchor 2 for C-terminal gap.")
        P = np.array([source_anchor1, source_anchor2])
        Q = np.array([target_anchor1, target_anchor2])
    else:
        raise ValueError("Unknown terminal type specified. Use 'N' or 'C'.")
        
    R_local, t_local = kabsch(P, Q)
    seg_atoms = []
    for key in segment_keys:
        seg_atoms.extend(folded_res_dict[key])
    adjusted_atoms = []
    for atom in seg_atoms:
        new_atom = atom.copy()
        coord = np.array([atom["x"], atom["y"], atom["z"]])
        new_coord = np.dot(R_local, coord) + t_local
        new_atom["x"], new_atom["y"], new_atom["z"] = new_coord.tolist()
        adjusted_atoms.append(new_atom)
    transformed_points = np.dot(R_local, P.T).T + t_local
    local_rmsd = np.sqrt(np.mean(np.sum((transformed_points - Q)**2, axis=1)))
    return adjusted_atoms, local_rmsd

########################################
# 5. Transfer Missing Segments (Entirely Missing Residues Only)
########################################

def oldtransfer_missing_segments(partial_atoms, folded_atoms):
    """
    For each chain, identify contiguous missing residue segments (residues present in folded but absent in partial).
    For each segment, perform local stitching using the appropriate routine (terminal or internal)
    and insert the transformed backbone atoms into the partial structure.
    Returns:
      - combined_atoms: list of atoms representing the merged structure.
      - transferred_events: list of dicts with chain, start, end, and local RMSD.
    """
    partial_res = group_by_residue(partial_atoms)
    folded_res = group_by_residue(folded_atoms)
    partial_keys = set(partial_res.keys())
    folded_keys = list(folded_res.keys())
    missing_segments = group_missing_segments(folded_keys, partial_keys)
    combined_atoms = list(partial_atoms)
    transferred_events = []
    for chain, segments in missing_segments.items():
        for segment in segments:
            if (chain, segment[0][1] - 1, "") not in partial_res:
                adjusted_atoms, local_rmsd = local_stitch_terminal(segment, folded_res, partial_res)
            else:
                adjusted_atoms, local_rmsd = local_stitch_internal(segment, folded_res, partial_res)
            event = {"chain": chain, "start": segment[0][1], "end": segment[-1][1], "rmsd": local_rmsd}
            transferred_events.append(event)
            group_adj = group_by_residue(adjusted_atoms)
            for key in segment:
                if key in group_adj:
                    bb = get_backbone_atoms(group_adj[key])
                    # For terminal segments, insert CA and N; for internal, insert CA for middle and CA+endpoint atom for endpoints.
                    if key == segment[0]:
                        if "CA" in bb:
                            new_atom = bb["CA"].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)
                        if "N" in bb:
                            new_atom = bb["N"].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)
                    elif key == segment[-1]:
                        if "CA" in bb:
                            new_atom = bb["CA"].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)
                        if "C" in bb:
                            new_atom = bb["C"].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)
                    else:
                        if "CA" in bb:
                            new_atom = bb["CA"].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)
    return combined_atoms, transferred_events


def transfer_missing_segments(partial_atoms, folded_atoms):
    """
    For each chain, identify contiguous missing residue segments (residues present in folded but absent in partial).
    For each segment, perform local stitching using the appropriate routine (terminal or internal)
    and insert the transformed backbone atoms into the partial structure.
    Returns:
      - combined_atoms: list of atoms representing the merged structure.
      - transferred_events: list of dicts with chain, start, end, and local RMSD.
    """
    partial_res = group_by_residue(partial_atoms)
    folded_res = group_by_residue(folded_atoms)
    partial_keys = set(partial_res.keys())
    folded_keys = list(folded_res.keys())
    missing_segments = group_missing_segments(folded_keys, partial_keys)
    combined_atoms = list(partial_atoms)
    transferred_events = []
    for chain, segments in missing_segments.items():
        for segment in segments:
            # Determine if the gap is terminal:
            if (chain, segment[0][1] - 1, "") not in partial_res:
                # Missing left anchor: N-terminal gap.
                adjusted_atoms, local_rmsd = local_stitch_terminal(segment, folded_res, partial_res, terminal='N')
            elif (chain, segment[-1][1] + 1, "") not in partial_res:
                # Missing right anchor: C-terminal gap.
                adjusted_atoms, local_rmsd = local_stitch_terminal(segment, folded_res, partial_res, terminal='C')
            else:
                adjusted_atoms, local_rmsd = local_stitch_internal(segment, folded_res, partial_res)
            event = {"chain": chain, "start": segment[0][1], "end": segment[-1][1], "rmsd": local_rmsd}
            transferred_events.append(event)
            group_adj = group_by_residue(adjusted_atoms)
            for key in segment:
                if key in group_adj:
                    bb = get_backbone_atoms(group_adj[key])
                    # For terminal segments, insert CA and N for N-terminal gaps and CA and C for C-terminal gaps;
                    # for internal segments, insert only CA for middle residues and CA+endpoint for endpoints.
                    if key == segment[0]:
                        if (chain, segment[0][1] - 1, "") not in partial_res:
                            # N-terminal gap: include CA and N.
                            if "CA" in bb:
                                new_atom = bb["CA"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                            if "N" in bb:
                                new_atom = bb["N"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                        elif (chain, segment[-1][1] + 1, "") not in partial_res:
                            # C-terminal gap: include N and CA.
                            if "N" in bb:
                                new_atom = bb["N"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                            if "CA" in bb:
                                new_atom = bb["CA"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                        else:
                            # Default: just add CA.
                            if "CA" in bb:
                                new_atom = bb["CA"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                    elif key == segment[-1]:
                        if (chain, segment[0][1] - 1, "") not in partial_res:
                            if "CA" in bb:
                                new_atom = bb["CA"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                            if "C" in bb:
                                new_atom = bb["C"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                        elif (chain, segment[-1][1] + 1, "") not in partial_res:
                            if "CA" in bb:
                                new_atom = bb["CA"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                            if "C" in bb:
                                new_atom = bb["C"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                        else:
                            if "CA" in bb:
                                new_atom = bb["CA"].copy()
                                new_atom["record_name"] = "ATOM"
                                combined_atoms.append(new_atom)
                    else:
                        if "CA" in bb:
                            new_atom = bb["CA"].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)
    return combined_atoms, transferred_events

########################################
# 6. Folding with ESM3 (using provided modules)
########################################

def fold_sequence(sequence, output_pdb):
    """
    Fold the given protein sequence using ESM3 and write the folded structure to output_pdb.
    Expects that ESM3 and torch are installed.
    """
    try:
        from esm.models.esm3 import ESM3
        from esm.sdk.api import ESMProtein, GenerationConfig
        from esm.pretrained import load_local_model
        from esm.utils.constants.models import normalize_model_name
        import torch
    except ImportError as e:
        print("Error: ESM3 modules not found. Ensure ESM3 and torch are installed.")
        sys.exit(1)
    model_name = "esm3-open"
    model_name = normalize_model_name(model_name)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = load_local_model(model_name, device=device)
    if device.type != "cpu":
        model = model.to(torch.bfloat16)
    protein = ESMProtein(sequence=sequence)
    protein = model.generate(protein, GenerationConfig(track="structure",num_steps=len(sequence),temperature=1e-6))
    protein.to_pdb(output_pdb)
    print(f"Folded structure written to {output_pdb}")

########################################
# 7. Amide Bond Sanity Check
########################################

def check_amide_bonds(combined_atoms, transferred_events):
    """
    For each inserted segment (identified by its starting residue),
    if there is a preceding residue in the same chain, compute the distance
    between the C atom of the preceding residue and the N atom of the inserted residue.
    Typical peptide bond lengths are ~1.2–1.5 Å.
    """
    merged_res = group_by_residue(combined_atoms)
    for event in transferred_events:
        chain = event["chain"]
        start = event["start"]
        key_inserted = (chain, start, "")
        prev_key = (chain, start - 1, "")
        if prev_key in merged_res and key_inserted in merged_res:
            prev_bb = get_backbone_atoms(merged_res[prev_key])
            curr_bb = get_backbone_atoms(merged_res[key_inserted])
            if "C" in prev_bb and "N" in curr_bb:
                C_coord = np.array([prev_bb["C"]["x"], prev_bb["C"]["y"], prev_bb["C"]["z"]])
                N_coord = np.array([curr_bb["N"]["x"], curr_bb["N"]["y"], curr_bb["N"]["z"]])
                bond_length = np.linalg.norm(C_coord - N_coord)
                print(f"Amide bond between {prev_key} and {key_inserted}: {bond_length:.3f} Å")
                if not (1.2 <= bond_length <= 1.5):
                    print(f"Warning: Unusual amide bond length for residue {key_inserted}: {bond_length:.3f} Å")
            else:
                print(f"Could not check amide bond for residue {key_inserted} (missing C or N).")

########################################
# 8. Main Workflow
########################################

def fold_and_transfer(partial_pdb, sequence, output_pdb):
    """
    Main workflow:
      1. Fold the sequence using ESM3 (write to a temporary folded PDB).
      2. Read the partial and folded PDB files.
      3. Compute the global alignment (using common CA atoms) and transform the folded structure.
      4. Report the global RMSD and write the globally fitted folded structure to "fitted_folded.pdb".
      5. Identify contiguous missing segments and perform local stitching using the chosen anchor atoms.
         Report the local RMSD for each segment (computed over the CA and O anchor correspondences).
      6. Insert the stitched segments (selected backbone atoms) into the partial structure.
      7. Perform an amide bond sanity check for the inserted segments.
      8. Write the final chimeric structure to the specified output file.
    """
    folded_temp = "folded_temp.pdb"
    # Step 1: Fold the sequence.
    fold_sequence(sequence, folded_temp)
    # Step 2: Read PDB files.
    try:
        partial_atoms = read_pdb(partial_pdb)
        folded_atoms = read_pdb(folded_temp)
        update_chain_id(folded_atoms,partial_atoms[0]["chain_id"])
    except Exception as e:
        print("Error reading PDB files:", e)
        sys.exit(1)

    # ★★★ DYNAMIC RENUMBERING OF THE FOLDED STRUCTURE ★★★
    # To ensure that common residues between partial and folded structures share the same keys,
    # we derive a chain-specific offset by mapping the partial structure’s sequence (from its residues)
    # onto the full sequence (provided as 'sequence').
    #
    # For example, suppose for chain A:
    #   - The partial structure’s first residue (by order) has experimental residue number R_exp.
    #   - Its one-letter code (from the atoms) appears in the full sequence at position i (1-indexed).
    # Then we define:
    #     offset = R_exp - i
    # and update each folded atom’s res_seq by adding offset.
    #
    # (This procedure works both when numbering is in the low single digits and when it is in the 300s.)
    partial_res = group_by_residue(partial_atoms)
    # For simplicity, assume a single chain; otherwise process each chain separately.
    chain = next(iter(partial_res))[0]
    partial_seq, partial_keys = extract_partial_sequence(partial_res, chain)
    indices = find_subsequence_indices(sequence, partial_seq)
    if indices is None or len(indices) == 0:
        print("Error: Could not align partial structure sequence with full sequence.")
        sys.exit(1)
    # The first residue in the partial structure (by order) corresponds to full sequence position (indices[0] + 1)
    print(partial_keys[0][1],end="->")
    print((indices[0] + 1))
    offset = partial_keys[0][1] - (indices[0] + 1)
    print(f"Computed renumbering offset for chain {chain}: {offset}")
    for atom in folded_atoms:
        atom["res_seq"] += offset
    # ★★★ end renumbering ★★★

    # Step 3: Global alignment using CA atoms.
    try:
        R, t, common_keys = align_folded_to_partial(partial_atoms, folded_atoms)
        print(f"Global alignment computed using {len(common_keys)} common residues.")
    except Exception as e:
        print("Error during global alignment:", e)
        sys.exit(1)
    transform_atoms(folded_atoms, R, t)
    rmsd_global = compute_alignment_rmsd(partial_atoms, folded_atoms, common_keys)
    print(f"Global alignment RMSD: {rmsd_global:.3f} Å")
    # Step 4: Write globally fitted folded structure.
    fitted_folded_file = "fitted_folded.pdb"
    write_pdb(fitted_folded_file, folded_atoms)
    print(f"Fitted folded structure written to {fitted_folded_file}")
    # Step 5 & 6: Identify missing segments and perform local stitching, then transfer.
    partial_res = group_by_residue(partial_atoms)
    folded_res = group_by_residue(folded_atoms)
    missing_segments = group_missing_segments(list(folded_res.keys()), set(partial_res.keys()))
    local_rmsd_reports = []
    combined_atoms, transferred_events = transfer_missing_segments(partial_atoms, folded_atoms)
    print(f"Inserted {len(transferred_events)} missing segments:")
    for event in transferred_events:
        print(f"  Chain {event['chain']} residues {event['start']}–{event['end']}, local RMSD = {event['rmsd']:.3f} Å")
        local_rmsd_reports.append(event)
    # Step 7: Amide bond check.
    check_amide_bonds(combined_atoms, transferred_events)
    # Step 8: Write final chimeric structure.
    try:
        write_pdb(output_pdb, combined_atoms)
        print(f"Final chimeric structure written to {output_pdb}")
    except Exception as e:
        print("Error writing output PDB:", e)
        sys.exit(1)
    try:
        os.remove(folded_temp)
    except Exception as e:
        print("Warning: could not remove temporary file:", folded_temp, e)
    if local_rmsd_reports:
        print("Local stitching RMSD reports:")
        for rep in local_rmsd_reports:
            print(f"  Chain {rep['chain']} residues {rep['start']}–{rep['end']}: RMSD = {rep['rmsd']:.3f} Å")
    else:
        print("No local stitching was performed.")

########################################
# 9. Command-Line Interface
########################################

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: fold_and_transfer.py <partial_pdb> <sequence> <output_pdb>")
        sys.exit(1)
    partial_pdb = sys.argv[1]
    sequence = sys.argv[2].strip()
    output_pdb = sys.argv[3]
    if not sequence:
        print("Error: Provided sequence is empty.")
        sys.exit(1)
    print("Starting folding, global alignment, local stitching, and segment transfer process...")
    fold_and_transfer(partial_pdb, sequence, output_pdb)

