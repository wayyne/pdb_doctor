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
"""

import numpy as np
import sys, os
from tqdm import tqdm
import warnings
import math

# Mute annoying warning
warnings.filterwarnings("ignore", message="Entity ID not found in metadata, using None as default")

# Store the original __init__ function before overriding
if not hasattr(tqdm, "_original_init"):
    tqdm._original_init = tqdm.__init__

# Define default tqdm settings globally
tqdm_defaults = {
    "ncols": 80,  # Fixed width
    "ascii": True,  # Use ASCII progress bars
    "leave": False   # Keep progress bar after completion
}

# Override tqdm globally with safe recursion handling
def custom_tqdm_init(self, *args, **kwargs):
    kwargs = {**tqdm_defaults, **kwargs}
    self._original_init(*args, **kwargs)

tqdm.__init__ = custom_tqdm_init

########################################
# Helper function for insertion code ordering
########################################

def insertion_order(i_code):
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
    keys = [k for k in res_dict if k[0] == chain]
    keys = sorted(keys, key=lambda k: (k[1], insertion_order(k[2])))
    seq = ""
    for key in keys:
        res_atoms = res_dict[key]
        res_name = res_atoms[0]["res_name"].upper()
        letter = three_to_one.get(res_name, 'X')
        seq += letter
    return seq, keys

def find_subsequence_indices(full_seq, subseq):
    import numpy as np
    len_full = len(full_seq)
    len_sub = len(subseq)

    dp = np.zeros((len_sub + 1, len_full + 1), dtype=int)
    backtrack = np.full((len_sub + 1, len_full + 1), -1, dtype=int)

    for i in range(1, len_sub + 1):
        for j in range(1, len_full + 1):
            if subseq[i - 1] == full_seq[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
                backtrack[i][j] = j - 1
            else:
                dp[i][j] = dp[i][j - 1]
                backtrack[i][j] = backtrack[i][j - 1]

    indices = []
    j = len_full
    for i in range(len_sub, 0, -1):
        j = backtrack[i][j]
        if j == -1:
            return None
        indices.append(j)
    return indices[::-1]

########################################
# 1. PDB File I/O Functions
########################################

def update_chain_id(atoms, new_chain_id):
    for atom in atoms:
        atom["chain_id"] = new_chain_id
    return atoms

def read_pdb(filename):
    atoms = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
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
                    "charge": charge
                }
                atoms.append(atom)
    return atoms

def write_pdb(filename, atoms):
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
    residues = {}
    for atom in atoms:
        if atom["record_name"] != "ATOM":
            continue
        key = (atom["chain_id"], atom["res_seq"], atom["i_code"])
        residues.setdefault(key, []).append(atom)
    return residues

def get_backbone_atoms(res_atoms):
    backbone = {}
    for atom in res_atoms:
        if atom["atom_name"] in ["N", "CA", "C", "O"]:
            backbone[atom["atom_name"]] = atom
    return backbone

def get_ca_coordinate(atom):
    return np.array([atom["x"], atom["y"], atom["z"]])

def get_atom_coordinate(atom):
    return np.array([atom["x"], atom["y"], atom["z"]])

########################################
# 3. Global Alignment (Kabsch using CA atoms)
########################################

def compute_centroid(coords):
    return np.mean(coords, axis=0)

def kabsch(P, Q):
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
    for atom in atoms:
        coord = np.array([atom["x"], atom["y"], atom["z"]])
        new_coord = np.dot(R, coord) + t
        atom["x"], atom["y"], atom["z"] = new_coord.tolist()

def compute_alignment_rmsd(partial_atoms, folded_atoms, common_keys):
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
        raise ValueError("Not enough common CA atoms for global alignment.")
    R, t = kabsch(np.array(folded_coords), np.array(partial_coords))
    return R, t, common_keys

########################################
# 4. Local Stitching Functions
########################################

def group_missing_segments(folded_keys, partial_keys):
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
                if (curr[1] - prev[1] == 1 and
                    insertion_order(curr[2]) - insertion_order(prev[2]) in [0, 1]):
                    current_seg.append(curr)
                else:
                    segs.append(current_seg)
                    current_seg = [curr]
            segs.append(current_seg)
            segments[chain] = segs
    return segments

BACKBONE_NAMES = ["N", "CA", "C", "O"]

def gather_anchor_atoms(res_dict_partial, res_dict_folded, key):
    partial_bb = get_backbone_atoms(res_dict_partial.get(key, []))
    folded_bb  = get_backbone_atoms(res_dict_folded.get(key, []))
    matched_pairs = []
    for name in BACKBONE_NAMES:
        if (name in partial_bb) and (name in folded_bb):
            matched_pairs.append((partial_bb[name], folded_bb[name]))
    return matched_pairs

def local_stitch_segment(segment_keys, folded_res_dict, partial_res_dict):
    if not segment_keys:
        raise ValueError("No residues in segment_keys.")
    chain = segment_keys[0][0]
    first_missing = segment_keys[0][1]
    last_missing  = segment_keys[-1][1]

    left_anchor_key  = (chain, first_missing - 1, "")
    right_anchor_key = (chain, last_missing + 1, "")

    anchor_pairs = []
    anchor_pairs += gather_anchor_atoms(partial_res_dict, folded_res_dict, left_anchor_key)
    anchor_pairs += gather_anchor_atoms(partial_res_dict, folded_res_dict, right_anchor_key)

    if len(anchor_pairs) < 2:
        raise ValueError("Not enough anchor-atom pairs to align the missing segment.")

    partial_coords = []
    folded_coords  = []
    for (p_atom, f_atom) in anchor_pairs:
        partial_coords.append([p_atom["x"], p_atom["y"], p_atom["z"]])
        folded_coords.append([f_atom["x"], f_atom["y"], f_atom["z"]])
    Q = np.array(partial_coords)
    P = np.array(folded_coords)
    R_local, t_local = kabsch(P, Q)

    seg_atoms = []
    for key in segment_keys:
        bb = get_backbone_atoms(folded_res_dict[key])
        seg_atoms.extend(bb.values())

    adjusted_atoms = []
    for atom in seg_atoms:
        new_atom = atom.copy()
        coord = np.array([atom["x"], atom["y"], atom["z"]])
        new_coord = R_local @ coord + t_local
        new_atom["x"], new_atom["y"], new_atom["z"] = new_coord.tolist()
        adjusted_atoms.append(new_atom)

    P_transformed = (R_local @ P.T).T + t_local
    local_rmsd = float(np.sqrt(np.mean(np.sum((P_transformed - Q)**2, axis=1))))
    return adjusted_atoms, local_rmsd

########################################
# 5. Transfer Missing Segments
########################################

def transfer_missing_segments(partial_atoms, folded_atoms):
    combined_atoms = list(partial_atoms)
    transferred_events = []

    partial_res = group_by_residue(partial_atoms)
    folded_res  = group_by_residue(folded_atoms)
    partial_keys = set(partial_res.keys())
    folded_keys  = list(folded_res.keys())

    missing_segments = group_missing_segments(folded_keys, partial_keys)

    for chain, segments in missing_segments.items():
        for segment in segments:
            partial_res = group_by_residue(combined_atoms)

            adjusted_atoms, local_rmsd = local_stitch_segment(segment, folded_res, partial_res)

            # Fix boundary bonds if needed 
            _fix_boundary_bonds(chain, segment, partial_res, adjusted_atoms)

            event = {"chain": chain, "start": segment[0][1], "end": segment[-1][1], "rmsd": local_rmsd}
            transferred_events.append(event)

            seg_dict = group_by_residue(adjusted_atoms)
            for key in segment:
                if key in seg_dict:
                    bb = get_backbone_atoms(seg_dict[key])
                    for name in ["N", "CA", "C", "O"]:
                        if name in bb:
                            new_atom = bb[name].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)

            partial_res = group_by_residue(combined_atoms)

    return combined_atoms, transferred_events

########################################
# 6. Folding with ESM3
########################################

def fold_sequence(sequence, output_pdb):
    try:
        from esm.models.esm3 import ESM3
        from esm.sdk.api import ESMProtein, GenerationConfig
        from esm.pretrained import load_local_model
        from esm.utils.constants.models import normalize_model_name
        import torch
    except ImportError:
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

def check_amide_bonds(combined_atoms, transferred_events, pdb_id):
    """
    For each inserted segment, report the left boundary AND the right boundary (if it exists).
    That is:
      - If there's a preceding residue (start-1), measure C->N
      - If there's a subsequent residue (end+1), measure C->N
    """
    merged_res = group_by_residue(combined_atoms)
    for event in transferred_events:
        chain = event["chain"]
        start = event["start"]
        end   = event["end"]

        # --- Left boundary check ---
        left_key = (chain, start - 1, "")
        inserted_key = (chain, start, "")
        if left_key in merged_res and inserted_key in merged_res:
            prev_bb = get_backbone_atoms(merged_res[left_key])
            curr_bb = get_backbone_atoms(merged_res[inserted_key])
            if "C" in prev_bb and "N" in curr_bb:
                C_coord = np.array([prev_bb["C"]["x"], prev_bb["C"]["y"], prev_bb["C"]["z"]])
                N_coord = np.array([curr_bb["N"]["x"], curr_bb["N"]["y"], curr_bb["N"]["z"]])
                bond_length = np.linalg.norm(C_coord - N_coord)
                print(f"Amide bond between {left_key} and {inserted_key}: {bond_length:.3f} Å")
                if not (1.2 <= bond_length <= 1.5):
                    print(f"Warning: Unusual amide bond length for PDBID {pdb_id} residue {inserted_key}: {bond_length:.3f} Å")
            else:
                print(f"Could not check left boundary bond for {inserted_key} (missing C or N).")

        # --- Right boundary check ---
        right_key = (chain, end + 1, "")
        inserted_end_key = (chain, end, "")
        if right_key in merged_res and inserted_end_key in merged_res:
            inserted_bb = get_backbone_atoms(merged_res[inserted_end_key])
            next_bb     = get_backbone_atoms(merged_res[right_key])
            if "C" in inserted_bb and "N" in next_bb:
                C_coord = np.array([inserted_bb["C"]["x"], inserted_bb["C"]["y"], inserted_bb["C"]["z"]])
                N_coord = np.array([next_bb["N"]["x"],  next_bb["N"]["y"],  next_bb["N"]["z"]])
                bond_length = np.linalg.norm(C_coord - N_coord)
                print(f"Amide bond between {inserted_end_key} and {right_key}: {bond_length:.3f} Å")
                if not (1.2 <= bond_length <= 1.5):
                    print(f"Warning: Unusual amide bond length for PDBID {pdb_id} residue {inserted_end_key}: {bond_length:.3f} Å")
            else:
                print(f"Could not check right boundary bond for {inserted_end_key} (missing C or N).")

########################################
# 8. Main Workflow
########################################

def fold_and_transfer(partial_pdb, sequence, output_pdb):
    folded_temp = "folded_temp.pdb"
    fold_sequence(sequence, folded_temp)

    try:
        partial_atoms = read_pdb(partial_pdb)
        folded_atoms  = read_pdb(folded_temp)
        update_chain_id(folded_atoms, partial_atoms[0]["chain_id"])
    except Exception as e:
        print("Error reading PDB files:", e)
        sys.exit(1)

    partial_res = group_by_residue(partial_atoms)
    chain = next(iter(partial_res))[0]
    partial_seq, partial_keys = extract_partial_sequence(partial_res, chain)
    indices = find_subsequence_indices(sequence, partial_seq)
    if not indices:
        print("Could not align partial structure with the full sequence.")
        sys.exit(1)
    offset = partial_keys[0][1] - (indices[0] + 1)
    for atom in folded_atoms:
        atom["res_seq"] += offset

    try:
        R, t, common_keys = align_folded_to_partial(partial_atoms, folded_atoms)
    except ValueError as e:
        print("Error in global alignment:", e)
        sys.exit(1)

    transform_atoms(folded_atoms, R, t)
    global_rmsd = compute_alignment_rmsd(partial_atoms, folded_atoms, common_keys)
    print(f"Global alignment RMSD: {global_rmsd:.3f} Å")

    fitted_folded_file = "fitted_folded.pdb"
    write_pdb(fitted_folded_file, folded_atoms)
    print(f"Fitted folded structure written to {fitted_folded_file}")

    combined_atoms, transferred_events = transfer_missing_segments(partial_atoms, folded_atoms)
    print(f"Inserted {len(transferred_events)} missing segments:")
    for evt in transferred_events:
        print(f"  Chain {evt['chain']} residues {evt['start']}–{evt['end']}, local RMSD = {evt['rmsd']:.3f} Å")

    check_amide_bonds(combined_atoms, transferred_events, partial_pdb[:4])

    write_pdb(output_pdb, combined_atoms)
    print(f"Final chimeric structure written to {output_pdb}")

    try:
        os.remove(folded_temp)
    except:
        pass

########################################
# 9. Boundary Fix Helpers
########################################

def _fix_boundary_bonds(chain, segment, partial_res, adjusted_atoms):
    """
    Adjusts peptide bonds at segment boundaries. If both boundaries are incorrect,
    it first corrects one boundary using translation, then rechecks and applies 
    rotation in a defined plane to correct the other boundary.
    """
    if not segment:
        return

    first_missing = segment[0][1]
    last_missing = segment[-1][1]

    left_key = (chain, first_missing - 1, "")
    right_key = (chain, last_missing + 1, "")

    seg_res = group_by_residue(adjusted_atoms)

    # Retrieve backbone atoms for left boundary
    partial_left_C = None
    new_left_N = None
    if left_key in partial_res:
        left_bb = get_backbone_atoms(partial_res[left_key])
        partial_left_C = left_bb.get("C", None)

    this_left = (chain, first_missing, "")
    if this_left in seg_res:
        left_bb_new = get_backbone_atoms(seg_res[this_left])
        new_left_N = left_bb_new.get("N", None)

    # Compute left bond distance
    left_dist = _bond_dist(partial_left_C, new_left_N) if (partial_left_C and new_left_N) else None

    # Retrieve backbone atoms for right boundary
    partial_right_N = None
    new_right_C = None
    if right_key in partial_res:
        right_bb = get_backbone_atoms(partial_res[right_key])
        partial_right_N = right_bb.get("N", None)

    this_right = (chain, last_missing, "")
    if this_right in seg_res:
        right_bb_new = get_backbone_atoms(seg_res[this_right])
        new_right_C = right_bb_new.get("C", None)

    # Compute right bond distance
    right_dist = _bond_dist(new_right_C, partial_right_N) if (partial_right_N and new_right_C) else None

    desired = 1.33  # Desired peptide bond length

    # Case 1: Only left boundary is incorrect → Translate left
    if left_dist and not right_dist:
        if not (1.2 <= left_dist <= 1.5):
            _translate_segment_for_boundary(partial_left_C, new_left_N, adjusted_atoms, desired)

    # Case 2: Only right boundary is incorrect → Translate right
    elif right_dist and not left_dist:
        if not (1.2 <= right_dist <= 1.5):
            _translate_segment_for_boundary(partial_right_N, new_right_C, adjusted_atoms, desired)

    # Case 3: Both boundaries are incorrect → Apply translation first, then rotation
    elif left_dist and right_dist:
        left_ok = (1.2 <= left_dist <= 1.5)
        right_ok = (1.2 <= right_dist <= 1.5)

        # Step 1: Apply translation to one boundary
        fixed_atom = None  # Atom that remains fixed for rotation later
        if not left_ok and not right_ok:
            # Prefer translating left boundary first
            _translate_segment_for_boundary(partial_left_C, new_left_N, adjusted_atoms, desired)
            left_ok = True
            fixed_atom = new_left_N  # The newly corrected atom

        # Step 2: Recompute the other boundary's bond distance
        if fixed_atom == new_left_N:
            right_dist = _bond_dist(new_right_C, partial_right_N) if (partial_right_N and new_right_C) else None
            right_ok = (1.2 <= right_dist <= 1.5)

        # Step 3: Apply rotation using the now corrected boundary
        if right_ok and not left_ok:
            _rotate_segment_in_plane(fixed_atom, partial_left_C, new_left_N, adjusted_atoms, desired)
        elif left_ok and not right_ok:
            _rotate_segment_in_plane(fixed_atom, partial_right_N, new_right_C, adjusted_atoms, desired)

def _bond_dist(a, b):
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return np.sqrt(dx*dx + dy*dy + dz*dz)

def _translate_segment_for_boundary(anchor, moved_atom, segment_atoms, desired_len):
    dist_now = _bond_dist(anchor, moved_atom)
    if dist_now < 1e-9:
        return
    scale = desired_len / dist_now
    dx = (moved_atom["x"] - anchor["x"]) * (scale - 1.0)
    dy = (moved_atom["y"] - anchor["y"]) * (scale - 1.0)
    dz = (moved_atom["z"] - anchor["z"]) * (scale - 1.0)
    for atm in segment_atoms:
        atm["x"] += dx
        atm["y"] += dy
        atm["z"] += dz

def _rotate_segment_in_plane(fixed_point, bond_atom1, bond_atom2, segment_atoms, desired_len):
    """
    Rotates the segment atoms around the fixed point in a plane defined by:
    - The two atoms forming the amide bond to be corrected.
    - The already corrected bond’s contributing atom (fixed_point).
    
    The function selects the optimal rotation angle to minimize bond length deviation.
    
    Parameters:
    - fixed_point: Dict with 'x', 'y', 'z' of the atom that remains fixed.
    - bond_atom1, bond_atom2: Dicts with 'x', 'y', 'z' for atoms forming the incorrect bond.
    - segment_atoms: List of atoms in the segment to be rotated.
    - desired_len: Target bond length for correction.
    """
    
    # Compute the normal vector to define the rotation plane
    vec1 = np.array([bond_atom1["x"] - fixed_point["x"],
                     bond_atom1["y"] - fixed_point["y"],
                     bond_atom1["z"] - fixed_point["z"]])
    
    vec2 = np.array([bond_atom2["x"] - fixed_point["x"],
                     bond_atom2["y"] - fixed_point["y"],
                     bond_atom2["z"] - fixed_point["z"]])
    
    plane_normal = np.cross(vec1, vec2)
    plane_normal /= np.linalg.norm(plane_normal)  # Normalize

    # Store original positions for restoring during testing
    original_coords = [(atm, (atm["x"], atm["y"], atm["z"])) for atm in segment_atoms]

    def restore_original():
        for (atm, (ox, oy, oz)) in original_coords:
            atm["x"], atm["y"], atm["z"] = ox, oy, oz

    # Optimization: Test multiple angles and select the best
    best_theta = None
    best_diff = float("inf")
    steps = 180  # Number of rotation steps

    for i in range(steps + 1):
        theta = 2 * math.pi * i / steps
        rotate_all_around_point(segment_atoms, fixed_point, plane_normal, theta)
        
        dtest = _bond_dist(bond_atom1, bond_atom2)
        diff = abs(dtest - desired_len)
        
        if diff < best_diff:
            best_diff = diff
            best_theta = theta
        
        restore_original()  # Reset for the next test

    # Apply best rotation if improvement is significant
    if best_theta is not None and best_diff < 0.3:
        rotate_all_around_point(segment_atoms, fixed_point, plane_normal, best_theta)


def _rotate_segment_around_line(fixA, fixB, offA, offB, segment_atoms, desired_len):
    import math
    v_line = np.array([fixB["x"] - fixA["x"],
                       fixB["y"] - fixA["y"],
                       fixB["z"] - fixA["z"]])
    axis_len = np.linalg.norm(v_line)
    if axis_len < 1e-8:
        return
    original_coords = [(atm, (atm["x"], atm["y"], atm["z"])) for atm in segment_atoms]

    def restore_original():
        for (atm, (ox, oy, oz)) in original_coords:
            atm["x"], atm["y"], atm["z"] = ox, oy, oz

    best_theta = None
    best_diff = 9999.0
    steps = 180
    for i in range(steps+1):
        theta = 2*math.pi*i/steps
        rotate_all_around_line(segment_atoms, fixA, fixB, theta)
        dtest = _bond_dist(offA, offB)
        diff = abs(dtest - desired_len)
        if diff < best_diff:
            best_diff = diff
            best_theta = theta
        restore_original()

    if best_theta is not None and best_diff < 0.3:
        rotate_all_around_line(segment_atoms, fixA, fixB, best_theta)

def rotate_all_around_point(segment_atoms, fixed_point, plane_normal, theta):
    """
    Rotates all atoms in `segment_atoms` around a `fixed_point` within a plane 
    defined by its normal vector `plane_normal` by an angle `theta` (radians).
    
    Parameters:
    - segment_atoms: List of atom dictionaries with 'x', 'y', 'z' coordinates.
    - fixed_point: Dict with 'x', 'y', 'z' defining the rotation center.
    - plane_normal: Normalized numpy array (3,) defining the plane of rotation.
    - theta: Rotation angle in radians.
    """

    # Ensure the normal is a unit vector
    plane_normal = plane_normal / np.linalg.norm(plane_normal)

    # Rodrigues' rotation formula components
    K = np.array([
        [0, -plane_normal[2], plane_normal[1]],
        [plane_normal[2], 0, -plane_normal[0]],
        [-plane_normal[1], plane_normal[0], 0]
    ])
    
    I = np.eye(3)
    R = I + math.sin(theta) * K + (1 - math.cos(theta)) * (K @ K)  # Rotation matrix

    # Apply rotation to each atom in the segment
    for atom in segment_atoms:
        # Convert to vector relative to the fixed point
        v = np.array([atom["x"] - fixed_point["x"],
                      atom["y"] - fixed_point["y"],
                      atom["z"] - fixed_point["z"]])

        # Apply rotation
        v_rot = R @ v

        # Update atom coordinates
        atom["x"] = fixed_point["x"] + v_rot[0]
        atom["y"] = fixed_point["y"] + v_rot[1]
        atom["z"] = fixed_point["z"] + v_rot[2]


def rotate_all_around_line(segment_atoms, p1, p2, theta):
    axis = np.array([p2["x"] - p1["x"],
                     p2["y"] - p1["y"],
                     p2["z"] - p1["z"]])
    norm_axis = np.linalg.norm(axis)
    if norm_axis < 1e-9:
        return
    axis /= norm_axis

    px, py, pz = p1["x"], p1["y"], p1["z"]
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    ux, uy, uz = axis

    for atm in segment_atoms:
        x0 = atm["x"] - px
        y0 = atm["y"] - py
        z0 = atm["z"] - pz
        dot = ux*x0 + uy*y0 + uz*z0
        crossx = uy*z0 - uz*y0
        crossy = uz*x0 - ux*z0
        crossz = ux*y0 - uy*x0

        rx = x0*cos_t + crossx*sin_t + ux*dot*(1 - cos_t)
        ry = y0*cos_t + crossy*sin_t + uy*dot*(1 - cos_t)
        rz = z0*cos_t + crossz*sin_t + uz*dot*(1 - cos_t)

        atm["x"] = rx + px
        atm["y"] = ry + py
        atm["z"] = rz + pz

########################################
# 10. Command-Line Interface
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

    print("=" * 80)
    print("Starting folding, global alignment, local stitching, and segment transfer process...")
    fold_and_transfer(partial_pdb, sequence, output_pdb)

