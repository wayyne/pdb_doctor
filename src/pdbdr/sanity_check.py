"""
Module: sanity_check
Author: [Your Name]
Date: [YYYY-MM-DD]
Description:
    This module provides a sanity check for stitched protein structures.
    In particular, it verifies that the amide bonds at the boundaries of
    transferred segments fall within the expected distance range.
"""

import numpy as np

from .pdb_io import group_by_residue, get_backbone_atoms
from .constants import VDW_RADII


def check_amide_bonds(combined_atoms, transferred_events, pdb_id):
    """
    Check the amide bond lengths at the boundaries of transferred segments.

    This function examines the left and right boundaries of each transferred
    segment by computing the distance between the C atom of the preceding
    residue and the N atom of the inserted residue (and vice versa for the
    right boundary). It prints the bond lengths and issues a warning if the
    bond length falls outside the acceptable range (1.2 Å to 1.5 Å).

    Args:
        combined_atoms (list): List of atom dictionaries representing the
            merged structure.
        transferred_events (list): List of dictionaries describing each
            transfer event. Each event must include:
                - "chain": Chain identifier.
                - "start": Starting residue number of the inserted segment.
                - "end": Ending residue number of the inserted segment.
        pdb_id (str): PDB identifier (used for reporting warnings).
    """
    # Group atoms by residue using a common key (chain, res_seq, i_code)
    merged_res = group_by_residue(combined_atoms)

    for event in transferred_events:
        chain = event["chain"]
        start = event["start"]
        end = event["end"]

        # --- Left Boundary Check ---
        left_key = (chain, start - 1, "")
        inserted_key = (chain, start, "")
        if left_key in merged_res and inserted_key in merged_res:
            prev_bb = get_backbone_atoms(merged_res[left_key])
            curr_bb = get_backbone_atoms(merged_res[inserted_key])
            if "C" in prev_bb and "N" in curr_bb:
                C_coord = np.array([
                    prev_bb["C"]["x"],
                    prev_bb["C"]["y"],
                    prev_bb["C"]["z"]
                ])
                N_coord = np.array([
                    curr_bb["N"]["x"],
                    curr_bb["N"]["y"],
                    curr_bb["N"]["z"]
                ])
                bond_length = np.linalg.norm(C_coord - N_coord)
                print(f"Amide bond between {left_key} and {inserted_key}: "
                      f"{bond_length:.3f} Å")
                if not (1.2 <= bond_length <= 1.5):
                    print(f"Warning: Unusual amide bond length for PDBID {pdb_id} "
                          f"residue {inserted_key}: {bond_length:.3f} Å")
            else:
                print(f"Could not check left boundary bond for {inserted_key} "
                      f"(missing C or N).")

        # --- Right Boundary Check ---
        right_key = (chain, end + 1, "")
        inserted_end_key = (chain, end, "")
        if right_key in merged_res and inserted_end_key in merged_res:
            inserted_bb = get_backbone_atoms(merged_res[inserted_end_key])
            next_bb = get_backbone_atoms(merged_res[right_key])
            if "C" in inserted_bb and "N" in next_bb:
                C_coord = np.array([
                    inserted_bb["C"]["x"],
                    inserted_bb["C"]["y"],
                    inserted_bb["C"]["z"]
                ])
                N_coord = np.array([
                    next_bb["N"]["x"],
                    next_bb["N"]["y"],
                    next_bb["N"]["z"]
                ])
                bond_length = np.linalg.norm(C_coord - N_coord)
                print(f"Amide bond between {inserted_end_key} and {right_key}: "
                      f"{bond_length:.3f} Å")
                if not (1.2 <= bond_length <= 1.5):
                    print(f"Warning: Unusual amide bond length for PDBID {pdb_id} "
                          f"residue {inserted_end_key}: {bond_length:.3f} Å")
            else:
                print(f"Could not check right boundary bond for "
                      f"{inserted_end_key} (missing C or N).")


def check_for_clashes(added_atoms, experimental_atoms, tolerance=0.0):
    """
    Check for steric clashes between newly added atoms (from the filled/stitched
    structure) and the experimental (anchor) atoms from the partial structure.
    
    This function is meant to be run after the filled structure has been aligned
    and its inserted atoms have been added to the partial structure. It uses
    vectorized computations to calculate pairwise distances between the inserted
    atoms and the experimental atoms. A clash is recorded if the distance between
    any inserted atom and any experimental atom is less than the sum of their
    van der Waals radii (minus an optional tolerance).
    
    Note:
        This function intentionally ignores any atoms that are part of the
        experimentally resolved (anchor) regions. It is assumed that the caller
        has separated the inserted atoms from the experimental ones.
    
    Args:
        added_atoms (list): List of atom dictionaries that were inserted from the
            folded/filled structure.
        experimental_atoms (list): List of atom dictionaries from the partial,
            experimentally determined structure.
        tolerance (float): Tolerance (in Å) to subtract from the sum of van der Waals
            radii (default is 0.0).
    
    Returns:
        list: A list of tuples (added_atom, experimental_atom, distance, overlap)
              for each detected clash, where:
                - added_atom: The inserted atom dictionary.
                - experimental_atom: The conflicting experimental atom.
                - distance: The Euclidean distance between them.
                - overlap: How much the two atoms overlap (in Å).
    """
    if not added_atoms or not experimental_atoms:
        return []

    # Convert positions to numpy arrays.
    added_positions = np.array([
        [atom["x"], atom["y"], atom["z"]] for atom in added_atoms
    ])
    exp_positions = np.array([
        [atom["x"], atom["y"], atom["z"]] for atom in experimental_atoms
    ])

    # Get van der Waals radii (defaulting to 1.7 Å if not found).
    default_radius = 1.7
    added_radii = np.array([
        VDW_RADII.get(atom["element"].upper(), default_radius) for atom in added_atoms
    ])
    exp_radii = np.array([
        VDW_RADII.get(atom["element"].upper(), default_radius) for atom in experimental_atoms
    ])

    # Compute pairwise Euclidean distances.
    diff = added_positions[:, None, :] - exp_positions[None, :, :]
    dists = np.linalg.norm(diff, axis=2)  # Shape: (n_added, n_exp)

    # Compute allowed minimum distances for each pair.
    allowed = added_radii[:, None] + exp_radii[None, :] - tolerance

    # Identify pairs where the distance is below the allowed value.
    clash_mask = dists < allowed
    i_new, i_exp = np.where(clash_mask)

    clashes = []
    for i, j in zip(i_new, i_exp):
        distance = dists[i, j]
        overlap = allowed[i, j] - distance
        clashes.append((added_atoms[i], experimental_atoms[j], distance, overlap))

    return clashes

def old_check_structure_for_clashes(combined_atoms, transferred_events, tolerance=0.0):
    """
    Check the entire final structure (combined_atoms) for steric clashes at the
    residue level. This function ignores clashes between covalently bonded
    residues (i.e. adjacent residues in the same chain) and then assigns any
    detected clash to a transferred segment (if either residue is within that
    segment's residue range).

    For each residue pair (not directly bonded) in the structure, the function
    computes all pairwise distances between atoms. A clash is flagged if any
    distance is less than the sum of the van der Waals radii (minus tolerance).

    Args:
        combined_atoms (list): List of atom dictionaries representing the final structure.
        transferred_events (list): List of dictionaries for transferred segments.
            Each dictionary must have keys "chain", "start", and "end".
        tolerance (float): Tolerance (in Å) to allow when comparing distances.
            Default is 0.3 Å.

    Returns:
        None. This function prints warnings for any transferred segment that has clashes,
        or an all-clear message if no clashes are detected.
    """
    from .pdb_io import group_by_residue

    # Group atoms by residue; each residue key is (chain, res_seq, i_code)
    res_groups = group_by_residue(combined_atoms)
    residue_list = []
    for res_key, atoms in res_groups.items():
        # Store a dictionary with the residue key and the positions of all atoms.
        positions = np.array([[atom["x"], atom["y"], atom["z"]] for atom in atoms])
        residue_list.append({
            "key": res_key,
            "positions": positions,
            "atoms": atoms
        })
    
    # Create a dictionary to collect clashes for each transferred segment.
    clashes_by_segment = { (ev["chain"], ev["start"], ev["end"]): [] 
                           for ev in transferred_events }
    
    n = len(residue_list)
    for i in range(n):
        for j in range(i + 1, n):
            res_i = residue_list[i]
            res_j = residue_list[j]
            key_i = res_i["key"]  # (chain, res_seq, i_code)
            key_j = res_j["key"]
            # If residues are in the same chain and directly consecutive,
            # skip the clash check (they are covalently bonded).
            if key_i[0] == key_j[0]:
                # If both insertion codes are empty and the difference is 1, skip.
                if key_i[2] == "" and key_j[2] == "" and abs(key_i[1] - key_j[1]) == 1:
                    continue

            # Compute all pairwise distances between atoms of these two residues.
            pos_i = res_i["positions"]  # shape (m, 3)
            pos_j = res_j["positions"]  # shape (n, 3)
            diff = pos_i[:, None, :] - pos_j[None, :, :]
            dists = np.linalg.norm(diff, axis=2)  # shape (m, n)
            
            # Build allowed distances for each atom pair.
            radii_i = np.array([VDW_RADII.get(atom["element"].upper(), 1.7) for atom in res_i["atoms"]])
            radii_j = np.array([VDW_RADII.get(atom["element"].upper(), 1.7) for atom in res_j["atoms"]])
            allowed = radii_i[:, None] + radii_j[None, :] - tolerance
            
            # If any pair has distance < allowed, record a clash.
            clash_mask = dists < allowed
            if not np.any(clash_mask):
                continue  # No clash between these two residues.

            # Compute minimal distance and average overlap for reporting.
            min_dist = np.min(dists[clash_mask])
            avg_overlap = np.mean(allowed[clash_mask] - dists[clash_mask])
            
            # Determine if either residue belongs to a transferred segment.
            # A residue is considered part of a transferred segment if its res_seq is between
            # the start and end of the segment (and on the same chain).
            for ev in transferred_events:
                if ev["chain"] != key_i[0]:
                    continue
                in_seg_i = (ev["start"] <= key_i[1] <= ev["end"])
                in_seg_j = (ev["start"] <= key_j[1] <= ev["end"])
                if in_seg_i or in_seg_j:
                    seg_key = (ev["chain"], ev["start"], ev["end"])
                    clashes_by_segment[seg_key].append((key_i, key_j, min_dist, avg_overlap))
    
    # Report the results.
    any_clash = False
    for seg_key, clashes in clashes_by_segment.items():
        if clashes:
            any_clash = True
            print(f"[WARNING] Detected {len(clashes)} clashes in transferred segment {seg_key}:")
            for clash in clashes:
                print(f"  Residues {clash[0]} and {clash[1]}: min distance {clash[2]:.2f} Å, avg overlap {clash[3]:.2f} Å")
    if not any_clash:
        print("All clear: no clashes detected in the combined structure.")



def check_structure_for_clashes(combined_atoms, transferred_events, tolerance=0.0):
    """
    Check the entire final structure for steric clashes at the residue level.
    
    For each pair of residues in the combined structure, the function iterates over
    every pair of heavy atoms (ignoring hydrogens). If the measured distance is less
    than the sum of the van der Waals radii (minus a tolerance), a clash is recorded.
    
    When the two residues are in the same chain and are directly consecutive, the
    expected covalent amide bond (between the "C" atom of the first and the "N" atom
    of the next) is ignored.
    
    Each detected clash is printed with full details: residue identifiers, atom names,
    measured distance, allowed distance, and the overlap amount. In addition, each clash
    is assigned to a transferred segment if either residue belongs to that segment.
    Clashes outside any transferred segment are reported as "global" clashes.
    
    Args:
        combined_atoms (list): List of atom dictionaries representing the final structure.
        transferred_events (list): List of transferred segment dictionaries. Each dictionary
            must include keys "chain", "start", and "end" (residue numbers).
        tolerance (float): Tolerance (in Å) to subtract from the sum of van der Waals radii.
            Default is 0.4 Å.
    
    Returns:
        None. Detailed clash information is printed.
    """
    # Group atoms by residue; each residue key is (chain, res_seq, i_code)
    res_groups = group_by_residue(combined_atoms)
    
    # Build a sorted list of residues.
    residues = []
    for key, atoms in res_groups.items():
        residues.append({"key": key, "atoms": atoms})
    residues.sort(key=lambda r: (r["key"][0], r["key"][1], r["key"][2]))
    
    # Prepare a dictionary to accumulate clashes per transferred segment.
    clashes_by_segment = { (ev["chain"], ev["start"], ev["end"]): [] for ev in transferred_events }
    global_key = ("global", 0, 0)
    clashes_by_segment[global_key] = []  # For clashes not belonging to any segment.
    
    # Compare each residue pair.
    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            res_i = residues[i]
            res_j = residues[j]
            key_i = res_i["key"]  # (chain, res_seq, i_code)
            key_j = res_j["key"]
            
            # If both residues are from the same chain and are consecutive (with empty insertion codes),
            # then ignore the expected covalent bond between residue i "C" and residue j "N".
            same_chain = (key_i[0] == key_j[0])
            adjacent = (key_i[2] == "" and key_j[2] == "" and abs(key_i[1] - key_j[1]) == 1)
            
            # We will still check other atom pairs in adjacent residues, but skip if the pair is exactly C vs N.
            # For residues not on the same chain or not consecutive, no special treatment.
            
            # Now check every pair of heavy atoms between these two residues.
            for atom_i in res_i["atoms"]:
                # Skip hydrogen atoms.
                if atom_i["element"].upper() == "H":
                    continue
                for atom_j in res_j["atoms"]:
                    if atom_j["element"].upper() == "H":
                        continue
                    # If the residues are adjacent in the same chain and this is the expected amide bond,
                    # skip the clash check.
                    if same_chain and adjacent:
                        if atom_i["atom_name"].strip() == "C" and atom_j["atom_name"].strip() == "N":
                            continue
                    pos_i = np.array([atom_i["x"], atom_i["y"], atom_i["z"]], dtype=float)
                    pos_j = np.array([atom_j["x"], atom_j["y"], atom_j["z"]], dtype=float)
                    dist = np.linalg.norm(pos_i - pos_j)
                    r_i = VDW_RADII.get(atom_i["element"].upper(), 1.7)
                    r_j = VDW_RADII.get(atom_j["element"].upper(), 1.7)
                    allowed = r_i + r_j - tolerance
                    if dist < allowed:
                        overlap = allowed - dist
                        clash_record = {
                            "residue_i": key_i,
                            "atom_i": atom_i["atom_name"].strip(),
                            "element_i": atom_i["element"].upper(),
                            "residue_j": key_j,
                            "atom_j": atom_j["atom_name"].strip(),
                            "element_j": atom_j["element"].upper(),
                            "distance": dist,
                            "allowed": allowed,
                            "overlap": overlap
                        }
                        # Determine if this clash belongs to a transferred segment.
                        assigned = False
                        for ev in transferred_events:
                            if ev["chain"] != key_i[0]:
                                continue
                            in_seg_i = (ev["start"] <= key_i[1] <= ev["end"])
                            in_seg_j = (ev["start"] <= key_j[1] <= ev["end"])
                            if in_seg_i or in_seg_j:
                                seg_key = (ev["chain"], ev["start"], ev["end"])
                                clashes_by_segment[seg_key].append(clash_record)
                                assigned = True
                        if not assigned:
                            clashes_by_segment[global_key].append(clash_record)
    
    # Report the clashes.
    any_clashes = False
    for seg_key, clashes in clashes_by_segment.items():
        if clashes:
            any_clashes = True
            if seg_key == global_key:
                print(f"[WARNING] Detected {len(clashes)} clashes in global regions (outside transferred segments):")
            else:
                print(f"[WARNING] Detected {len(clashes)} clashes in transferred segment {seg_key}:")
            for clash in clashes:
                print(f"  Residue {clash['residue_i']} ({clash['atom_i']}, {clash['element_i']}) and "
                      f"Residue {clash['residue_j']} ({clash['atom_j']}, {clash['element_j']}): "
                      f"distance = {clash['distance']:.2f} Å, allowed = {clash['allowed']:.2f} Å, "
                      f"overlap = {clash['overlap']:.2f} Å")
    if not any_clashes:
        print("All clear: no clashes detected in the combined structure.")

