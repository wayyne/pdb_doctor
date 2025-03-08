"""
Module: stitching
Author: [Your Name]
Date: [YYYY-MM-DD]
Description:
    This module performs stitching of partial protein structure segments
    with filled regions. It provides functions to identify missing segments,
    gather anchor residues, fix boundary bonds, bridge missing residues,
    check for steric clashes, and transfer missing segments into a
    partial structure.
"""

import numpy as np
import math
from copy import deepcopy

from .pdb_io import (group_by_residue, get_backbone_atoms, insertion_order)
from .alignment import kabsch
from .sanity_check import check_for_clashes
from .constants import (
    BACKBONE_NAMES,
    BOND_N_CA, BOND_CA_C, BOND_C_O,
    ANGLE_NCAC, ANGLE_CCO, THREE_TO_ONE
)


def group_missing_segments(filled_keys, partial_keys):
    """
    Group consecutive missing residue keys into segments.

    Args:
        filled_keys (iterable): Residue keys from the filled structure.
            Each key is expected to be a tuple (chain, res_seq, i_code).
        partial_keys (set): Residue keys present in the partial structure.

    Returns:
        dict: Mapping from chain identifier to a list of segments. Each
            segment is a list of residue keys that are missing.
    """
    segments = {}
    chain_keys = {}
    # Group keys by chain.
    for key in filled_keys:
        chain_keys.setdefault(key[0], []).append(key)
    for chain, keys in chain_keys.items():
        # Sort keys by residue sequence and insertion code.
        keys = sorted(keys, key=lambda k: (k[1], insertion_order(k[2])))
        # Identify keys missing from the partial structure.
        missing = [k for k in keys if k not in partial_keys]
        segs = []
        if missing:
            current_seg = [missing[0]]
            # Group consecutive missing keys into segments.
            for prev, curr in zip(missing, missing[1:]):
                if (curr[1] - prev[1] == 1 and
                        insertion_order(curr[2]) - insertion_order(prev[2])
                        in [0, 1]):
                    current_seg.append(curr)
                else:
                    segs.append(current_seg)
                    current_seg = [curr]
            segs.append(current_seg)
            segments[chain] = segs
    return segments


def gather_anchor_residues(res_dict_partial, res_dict_filled, chain,
                           start_res_seq, window_size):
    """
    Gather matching backbone atom pairs between partial and filled
    residue dictionaries for anchoring.

    Args:
        res_dict_partial (dict): Mapping of residue keys to atom lists from
            the partial structure.
        res_dict_filled (dict): Mapping of residue keys to atom lists from
            the filled structure.
        chain (str): Chain identifier.
        start_res_seq (int): Starting residue sequence number.
        window_size (int): Number of consecutive residues to consider.

    Returns:
        list: A list of tuples, each containing a pair of backbone atom
            dictionaries (partial, filled) for matching residue positions.
    """
    matched_pairs = []
    for offset in range(window_size + 1):
        candidate_key = (chain, start_res_seq + offset, "")
        if (candidate_key in res_dict_partial and
                candidate_key in res_dict_filled):
            partial_bb = get_backbone_atoms(res_dict_partial[candidate_key])
            filled_bb = get_backbone_atoms(res_dict_filled[candidate_key])
            for name in ["N", "CA", "C", "O"]:
                if name in partial_bb and name in filled_bb:
                    matched_pairs.append((partial_bb[name], filled_bb[name]))
    return matched_pairs


def check_boundary_amide(res_left_atoms, res_right_atoms):
    """
    Check if the amide bond between two residues is within acceptable
    bounds.

    Args:
        res_left_atoms (list or dict): Atom information for the left residue.
        res_right_atoms (list or dict): Atom information for the right residue.

    Returns:
        tuple: (bool, float or None) where the bool indicates whether the bond
            length is acceptable (typically 1.2-1.5 Å) and the float is the
            computed distance (or None if not available).
    """
    left_bb = get_backbone_atoms(res_left_atoms)
    right_bb = get_backbone_atoms(res_right_atoms)

    if "C" not in left_bb or "N" not in right_bb:
        return True, None

    c_coord = np.array([left_bb["C"]["x"],
                        left_bb["C"]["y"],
                        left_bb["C"]["z"]])
    n_coord = np.array([right_bb["N"]["x"],
                        right_bb["N"]["y"],
                        right_bb["N"]["z"]])
    dist = np.linalg.norm(c_coord - n_coord)
    if 1.2 <= dist <= 1.5:
        return True, dist
    else:
        return False, dist


def _bond_dist(a, b):
    """
    Compute the Euclidean distance between two atoms.

    Args:
        a (dict): Atom with keys 'x', 'y', 'z'.
        b (dict): Atom with keys 'x', 'y', 'z'.

    Returns:
        float or None: The distance between a and b, or None if either is
            None.
    """
    if a is None or b is None:
        return None
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return np.sqrt(dx * dx + dy * dy + dz * dz)


def _translate_segment_for_boundary(anchor, moved_atom, segment_atoms,
                                    desired_len):
    """
    Translate a group of atoms so that the bond length between the anchor
    and the moved_atom becomes the desired length.

    Args:
        anchor (dict): Reference atom.
        moved_atom (dict): Atom to be adjusted.
        segment_atoms (list): List of atom dictionaries to translate.
        desired_len (float): Target bond length.
    """
    dist_now = _bond_dist(anchor, moved_atom)
    if dist_now is None or dist_now < 1e-9:
        return
    scale = desired_len / dist_now
    dx = (moved_atom["x"] - anchor["x"]) * (scale - 1.0)
    dy = (moved_atom["y"] - anchor["y"]) * (scale - 1.0)
    dz = (moved_atom["z"] - anchor["z"]) * (scale - 1.0)
    for atm in segment_atoms:
        atm["x"] += dx
        atm["y"] += dy
        atm["z"] += dz


def _fix_boundary_bonds(chain, segment, partial_res, adjusted_atoms):
    """
    Adjust the boundary bonds of a stitched segment so that the N-C amide
    bonds match the desired distances.

    Args:
        chain (str): Chain identifier.
        segment (list): List of residue keys (tuples) in the segment.
        partial_res (dict): Partial structure residue mapping.
        adjusted_atoms (list): List of adjusted atom dictionaries (to be modified).
    """
    if not segment:
        return
    first_missing = segment[0][1]
    last_missing = segment[-1][1]

    left_key = (chain, first_missing - 1, "")
    right_key = (chain, last_missing + 1, "")

    seg_res = group_by_residue(adjusted_atoms)

    partial_left_C = None
    new_left_N = None
    if left_key in partial_res:
        left_bb = get_backbone_atoms(partial_res[left_key])
        partial_left_C = left_bb.get("C", None)

    this_left = (chain, first_missing, "")
    if this_left in seg_res:
        left_bb_new = get_backbone_atoms(seg_res[this_left])
        new_left_N = left_bb_new.get("N", None)

    left_dist = _bond_dist(partial_left_C, new_left_N)

    partial_right_N = None
    new_right_C = None
    if right_key in partial_res:
        right_bb = get_backbone_atoms(partial_res[right_key])
        partial_right_N = right_bb.get("N", None)

    this_right = (chain, last_missing, "")
    if this_right in seg_res:
        right_bb_new = get_backbone_atoms(seg_res[this_right])
        new_right_C = right_bb_new.get("C", None)

    right_dist = _bond_dist(new_right_C, partial_right_N)

    desired = 1.33
    if left_dist and (not right_dist):
        if not (1.2 <= left_dist <= 1.5):
            _translate_segment_for_boundary(partial_left_C, new_left_N,
                                            adjusted_atoms, desired)
    elif right_dist and (not left_dist):
        if not (1.2 <= right_dist <= 1.5):
            _translate_segment_for_boundary(partial_right_N, new_right_C,
                                            adjusted_atoms, desired)
    elif left_dist and right_dist:
        # Both bonds exist; additional logic could be added here.
        pass


def _build_single_residue_between_anchors(left_anchor, right_anchor,
                                          desired_bond=1.33, scale=True):
    """
    Build a single residue that bridges two anchor atoms using ideal
    bond lengths and angles.

    Args:
        left_anchor (dict): Atom dictionary for the left anchor.
        right_anchor (dict): Atom dictionary for the right anchor.
        desired_bond (float): Desired bond length (default 1.33 Å).
        scale (bool): Whether to scale the ideal residue to match the gap.

    Returns:
        dict: A dictionary with keys 'N', 'CA', 'C', and 'O' representing the
            built residue atoms.
    """
    # Ideal bond lengths and angles (in Å and radians)
    bond_N_CA = 1.46
    bond_CA_C = 1.52
    bond_C_O = 1.23
    angle_NCAC = math.radians(111.0)
    angle_CCO = math.radians(120.0)

    # Define local coordinates for an ideal residue.
    N_local = np.array([0.0, 0.0, 0.0])
    CA_local = np.array([bond_N_CA, 0.0, 0.0])
    dx = bond_CA_C * math.cos(angle_NCAC)
    dy = bond_CA_C * math.sin(angle_NCAC)
    C_local = CA_local + np.array([dx, dy, 0.0])

    # Compute O atom position relative to C.
    v_C_CA = CA_local - C_local
    len_C_CA = np.linalg.norm(v_C_CA)
    if len_C_CA < 1e-9:
        O_local = C_local + np.array([0.0, -bond_C_O, 0.0])
    else:
        cosr = math.cos(angle_CCO)
        sinr = math.sin(angle_CCO)

        def rotate2d(vec, angle):
            x, y = vec[0], vec[1]
            xr = x * cosr - y * sinr
            yr = x * sinr + y * cosr
            return np.array([xr, yr, 0.0])

        v_c_o_2d = rotate2d(v_C_CA, angle_CCO)
        v_c_o_2d = v_c_o_2d / np.linalg.norm(v_c_o_2d) * bond_C_O
        O_local = C_local + v_c_o_2d

    local_NC = np.linalg.norm(C_local - N_local)
    L = np.array([left_anchor["x"],
                  left_anchor["y"],
                  left_anchor["z"]])
    R = np.array([right_anchor["x"],
                  right_anchor["y"],
                  right_anchor["z"]])
    anchor_dist = np.linalg.norm(R - L)
    needed_NC = anchor_dist - 2.0 * desired_bond
    if needed_NC < 0.5:
        needed_NC = 1.0
    factor = 1.0
    if scale and local_NC > 1e-9:
        factor = needed_NC / local_NC

    # Scale the local coordinates.
    N_loc_s = N_local * factor
    CA_loc_s = CA_local * factor
    C_loc_s = C_local * factor
    O_loc_s = O_local * factor

    # Align the ideal residue between anchors.
    direction_LR = R - L
    dist_LR = np.linalg.norm(direction_LR)
    if dist_LR < 1e-9:
        direction_LR = np.array([1.0, 0.0, 0.0])
    else:
        direction_LR /= dist_LR
    finalN = L + direction_LR * (desired_bond)
    finalC = R - direction_LR * (desired_bond)

    # Use Kabsch alignment on the N and C atoms.
    pLocal = np.vstack([N_loc_s, C_loc_s])
    pTarget = np.vstack([finalN, finalC])
    R_, t_ = kabsch(pLocal, pTarget)

    # Build the atom dictionaries for the residue.
    all_atoms = {
        "N":  deepcopy({"atom_name": "N",
                        "x": N_loc_s[0],
                        "y": N_loc_s[1],
                        "z": N_loc_s[2]}),
        "CA": deepcopy({"atom_name": "CA",
                        "x": CA_loc_s[0],
                        "y": CA_loc_s[1],
                        "z": CA_loc_s[2]}),
        "C":  deepcopy({"atom_name": "C",
                        "x": C_loc_s[0],
                        "y": C_loc_s[1],
                        "z": C_loc_s[2]}),
        "O":  deepcopy({"atom_name": "O",
                        "x": O_loc_s[0],
                        "y": O_loc_s[1],
                        "z": O_loc_s[2]})
    }
    # Transform all atoms using the computed rotation and translation.
    for name, at in all_atoms.items():
        xyz = np.array([at["x"], at["y"], at["z"]])
        xyzT = R_.dot(xyz) + t_
        at["x"], at["y"], at["z"] = xyzT.tolist()
    return all_atoms


def _bridge_single_residue(res_atoms, left_anchor, right_anchor,
                           desired_bond=1.35, bond_tolerance=0.15,
                           max_iter=100):
    """
    Bridge a single residue between two anchors by iteratively adjusting
    atom positions to achieve the desired bond lengths.

    Args:
        res_atoms (dict): Dictionary of residue atoms.
        left_anchor (dict): Left anchor atom.
        right_anchor (dict): Right anchor atom.
        desired_bond (float): Target bond length (default 1.35 Å).
        bond_tolerance (float): Acceptable deviation from desired_bond.
        max_iter (int): Maximum number of iterations.

    Returns:
        list: A list of adjusted atom dictionaries.
    """
    # Use deepcopy to avoid modifying the original atoms.
    residue_atoms = [deepcopy(a) for a in res_atoms.values()]

    def get_coord(atom):
        return np.array([atom["x"], atom["y"], atom["z"]], dtype=float)

    def set_coord(atom, coord):
        atom["x"], atom["y"], atom["z"] = coord

    L_coord = get_coord(left_anchor)
    R_coord = get_coord(right_anchor)

    def dist(atom, anchor_coord):
        if atom is None:
            return None
        pos = get_coord(atom)
        return np.linalg.norm(pos - anchor_coord)

    def fix_left_bond(atoms, anchor_coord, nN_atom, desired):
        oldPos = get_coord(nN_atom)
        dist_now = np.linalg.norm(oldPos - anchor_coord)
        if dist_now < 1e-9:
            return
        shift_needed = desired - dist_now
        direction = (oldPos - anchor_coord) / dist_now
        shift_vec = shift_needed * direction
        for at in atoms:
            coord = get_coord(at)
            set_coord(at, coord + shift_vec)

    def fix_right_bond(atoms, anchor_coord, nC_atom, pivot_atom, desired):
        pivot = get_coord(pivot_atom)
        cpos = get_coord(nC_atom)
        dist_now = np.linalg.norm(cpos - anchor_coord)
        if dist_now < 1e-9:
            return
        v_pc = cpos - pivot
        v_pa = anchor_coord - pivot
        axis = np.cross(v_pc, v_pa)
        ax_len = np.linalg.norm(axis)
        if ax_len < 1e-9:
            fix_left_bond(atoms, anchor_coord, nC_atom, desired)
            return
        axis /= ax_len
        current_error = dist_now - desired
        if abs(current_error) < 1e-3:
            return
        step_deg = 2.0
        best_theta = 0
        best_dist = dist_now
        sign_options = [1, -1]

        def rodrigues_rotate(vec, axis, theta):
            c = math.cos(theta)
            s = math.sin(theta)
            crossA = np.cross(axis, vec)
            dotA = np.dot(axis, vec)
            return vec * c + crossA * s + axis * dotA * (1 - c)

        for sign in sign_options:
            for i in range(1, 6):
                theta = math.radians(step_deg * i * sign)
                shifted_cpos = cpos - pivot
                test_cpos = rodrigues_rotate(shifted_cpos, axis, theta) + pivot
                test_dist = np.linalg.norm(test_cpos - anchor_coord)
                if abs(test_dist - desired) < abs(best_dist - desired):
                    best_theta = theta
                    best_dist = test_dist
        if abs(best_theta) > 1e-9:
            c = math.cos(best_theta)
            s = math.sin(best_theta)
            crossA = np.array([[0, -axis[2], axis[1]],
                               [axis[2], 0, -axis[0]],
                               [-axis[1], axis[0], 0]], dtype=float)
            R_mat = np.eye(3) + s * crossA + (1 - c) * (crossA @ crossA)
            for at in atoms:
                pold = get_coord(at) - pivot
                pnew = R_mat @ pold
                set_coord(at, pivot + pnew)

    desired_min = desired_bond - bond_tolerance
    desired_max = desired_bond + bond_tolerance

    for _ in range(max_iter):
        nN_atom, nC_atom = None, None
        for a in residue_atoms:
            if a["atom_name"] == "N":
                nN_atom = a
            elif a["atom_name"] == "C":
                nC_atom = a
        if nN_atom is None or nC_atom is None:
            break
        fix_left_bond(residue_atoms, L_coord, nN_atom, desired_bond)
        fix_right_bond(residue_atoms, R_coord, nC_atom, nN_atom,
                       desired_bond)
        left_dist = dist(nN_atom, L_coord)
        right_dist = dist(nC_atom, R_coord)
        if (left_dist is not None and right_dist is not None and
                desired_min <= left_dist <= desired_max and
                desired_min <= right_dist <= desired_max):
            break
    return residue_atoms

def local_stitch_segment_extend(segment_keys, filled_res_dict, partial_res_dict, max_extend=5):
    """
    Stitch a missing segment by dynamically extending the missing region
    on either boundary if that boundary's amide bond is misaligned.
    
    Returns:
        (adjusted_atoms, local_rmsd, anchor_atoms, final_first, final_last, left_expanded, right_expanded)
        
      - adjusted_atoms: newly aligned backbone atoms for the final missing range
      - local_rmsd: RMSD of the final alignment
      - anchor_atoms: tuple (left_anchor_atoms, right_anchor_atoms) for reporting
      - final_first, final_last: the final expanded residue indices
      - left_expanded, right_expanded: how many times we extended on each boundary
    """
    import math
    import numpy as np
    from copy import deepcopy

    if not segment_keys:
        raise ValueError("segment_keys is empty.")

    chain = segment_keys[0][0]
    first_missing = segment_keys[0][1]
    last_missing = segment_keys[-1][1]

    # For final reporting, store the original anchor atoms
    left_anchor_key = (chain, first_missing - 1, "")
    right_anchor_key = (chain, last_missing + 1, "")
    left_anchor_atoms = partial_res_dict.get(left_anchor_key, None)
    right_anchor_atoms = partial_res_dict.get(right_anchor_key, None)
    anchor_atoms = (left_anchor_atoms, right_anchor_atoms)

    # Count how many anchors we have
    has_left_anchor = (left_anchor_key in partial_res_dict)
    has_right_anchor = (right_anchor_key in partial_res_dict)
    anchor_count = int(has_left_anchor) + int(has_right_anchor)

    print(f"\n=== local_stitch_segment_extend for {chain} {first_missing}-{last_missing} ===")
    print(f"  => anchor_count={anchor_count}")

    def build_segment_keys(start_res, end_res):
        return [(chain, r, "") for r in range(start_res, end_res + 1)]

    def gather_segment_atoms(skeys):
        out = []
        for sk in skeys:
            bb = get_backbone_atoms(filled_res_dict[sk])
            out.extend(list(bb.values()))
        return out

    def do_kabsch(seg_keys, left_k, right_k):
        """Return (adjusted_atoms, local_rmsd, left_dist, right_dist)."""
        anchor_pairs = []
        for akey in [left_k, right_k]:
            if akey in partial_res_dict and akey in filled_res_dict:
                p_bb = get_backbone_atoms(partial_res_dict[akey])
                f_bb = get_backbone_atoms(filled_res_dict[akey])
                for nm in ["N", "CA", "C", "O"]:
                    if nm in p_bb and nm in f_bb:
                        anchor_pairs.append((p_bb[nm], f_bb[nm]))

        seg_atoms_filled = gather_segment_atoms(seg_keys)
        if len(anchor_pairs) < 2:
            return [a.copy() for a in seg_atoms_filled], 999.9, None, None

        partial_coords = []
        filled_coords = []
        for (p_atom, f_atom) in anchor_pairs:
            partial_coords.append(np.array([p_atom["x"], p_atom["y"], p_atom["z"]], dtype=float))
            filled_coords.append(np.array([f_atom["x"], f_atom["y"], f_atom["z"]], dtype=float))

        P = np.array(filled_coords)   # "from"
        Q = np.array(partial_coords)  # "to"
        R_local, t_local = kabsch(P, Q)

        adjusted = []
        for atm in seg_atoms_filled:
            new_a = atm.copy()
            c = np.array([atm["x"], atm["y"], atm["z"]])
            new_c = R_local.dot(c) + t_local
            new_a["x"], new_a["y"], new_a["z"] = new_c.tolist()
            adjusted.append(new_a)

        # RMSD
        P_t = (R_local.dot(P.T)).T + t_local
        local_rmsd = float(math.sqrt(np.mean(np.sum((P_t - Q)**2, axis=1))))

        # boundary distances
        seg_res = group_by_residue(adjusted)
        left_dist = None
        right_dist = None

        if seg_keys and left_k in partial_res_dict:
            first_res = seg_keys[0]
            if first_res in seg_res:
                okL, distL = check_boundary_amide(
                    partial_res_dict[left_k],
                    seg_res[first_res]
                )
                left_dist = distL

        if seg_keys and right_k in partial_res_dict:
            last_res = seg_keys[-1]
            if last_res in seg_res:
                okR, distR = check_boundary_amide(
                    seg_res[last_res],
                    partial_res_dict[right_k]
                )
                right_dist = distR

        return adjusted, local_rmsd, left_dist, right_dist

    # single-res + both anchors => bridging
    if len(segment_keys) == 1 and anchor_count == 2:
        print("  => Single-res bridging with both anchors => _bridge_single_residue")
        single_key = segment_keys[0]
        filled_bb = get_backbone_atoms(filled_res_dict[single_key])
        left_bb = get_backbone_atoms(partial_res_dict[left_anchor_key])
        right_bb = get_backbone_atoms(partial_res_dict[right_anchor_key])
        if ("C" in left_bb and "N" in right_bb and
            "N" in filled_bb and "C" in filled_bb):
            bridged_atoms = _bridge_single_residue(filled_bb, left_bb["C"], right_bb["N"], 1.33)
            return bridged_atoms, 0.0, anchor_atoms, first_missing, last_missing, 0, 0

    # no anchors => raw copy
    if anchor_count == 0:
        print("  => No anchors => copy raw segment atoms.")
        seg_atoms = []
        for k in segment_keys:
            bb = get_backbone_atoms(filled_res_dict[k])
            seg_atoms.extend(list(bb.values()))
        return [a.copy() for a in seg_atoms], 0.0, anchor_atoms, first_missing, last_missing, 0, 0

    # single anchor => do single-boundary fix
    if anchor_count == 1:
        print("  => single-anchor => minimal alignment + _fix_boundary_bonds")
        seg_atoms = []
        for k in segment_keys:
            bb = get_backbone_atoms(filled_res_dict[k])
            seg_atoms.extend(list(bb.values()))
        adjusted = [a.copy() for a in seg_atoms]
        _fix_boundary_bonds(chain, segment_keys, partial_res_dict, adjusted)
        return adjusted, 0.0, anchor_atoms, first_missing, last_missing, 0, 0

    # anchor_count == 2 => internal => expansions
    current_first = first_missing
    current_last = last_missing

    best_atoms = []
    best_rmsd = 999.9
    done = False
    iteration = 0

    left_expanded = 0
    right_expanded = 0

    while iteration < 2 * max_extend + 10:
        iteration += 1
        seg_keys = build_segment_keys(current_first, current_last)
        left_k = (chain, current_first - 1, "")
        right_k = (chain, current_last + 1, "")

        print(f"\nIteration {iteration}: missing={current_first}-{current_last}, anchors={left_k}, {right_k}")
        (adj_atoms,
         local_rmsd,
         left_dist,
         right_dist) = do_kabsch(seg_keys, left_k, right_k)
        print(f"  => RMSD={local_rmsd:.3f}, left_dist={left_dist}, right_dist={right_dist}")

        if local_rmsd < best_rmsd:
            best_rmsd = local_rmsd
            best_atoms = adj_atoms

        # Check boundary: 1.2..1.5 is good
        left_good = True
        right_good = True
        if left_dist is not None:
            left_good = (1.2 <= left_dist <= 1.5)
        if right_dist is not None:
            right_good = (1.2 <= right_dist <= 1.5)

        print(f"  => left_good={left_good}, right_good={right_good}")

        if (left_dist is None and right_dist is None):
            # no valid anchor pairs
            print("  => no valid anchor pairs => stopping.")
            break

        # if both boundaries good => done
        if left_good and right_good:
            print("  => Both boundaries good => done")
            best_atoms = adj_atoms
            best_rmsd = local_rmsd
            done = True
            break

        # expand if boundary is bad
        expanded_this_round = False
        if not left_good and left_expanded < max_extend and (left_k in partial_res_dict):
            print(f"  => Incorporate left anchor {left_k} => expand missing range left")
            current_first -= 1
            left_expanded += 1
            expanded_this_round = True

        if not right_good and right_expanded < max_extend and (right_k in partial_res_dict):
            print(f"  => Incorporate right anchor {right_k} => expand missing range right")
            current_last += 1
            right_expanded += 1
            expanded_this_round = True

        # if we didn't expand => can't fix further
        if not expanded_this_round:
            print("  => can't expand further, stopping.")
            break

    if not done:
        print(f"  => final boundaries not perfect, best_rmsd={best_rmsd:.3f}")

    return best_atoms, best_rmsd, anchor_atoms, current_first, current_last, left_expanded, right_expanded


def great_local_stitch_segment(
    segment_keys,
    filled_res_dict,
    partial_res_dict,
    max_extend=5
):
    """
    Stitch a missing segment by dynamically extending the missing region
    on either boundary if that boundary is misaligned.

    How it works:
      1) Gather a "missing range" from first_missing..last_missing, with
         anchors (first_missing-1, last_missing+1) if they exist in partial.
      2) Kabsch align to these anchor keys (if 2 anchors). If 1 anchor, we do
         single-boundary logic. If 0 anchors, copy raw.
      3) Check left/right boundary bond lengths:
         - If left boundary is out of range, incorporate the left anchor
           residue into the missing segment. That means we shift the
           missing segment's start by -1, and the new left anchor becomes
           the old "anchor-1".
         - If right boundary is out of range, do the same on the right side,
           i.e. missing segment is extended by +1, and anchor becomes +2.
      4) Repeat the process, each time re-building the missing segment atoms
         from 'filled_res_dict' for the expanded range, until:
           - Both boundaries are acceptable, or
           - We can no longer extend (ran out of partial anchor residues),
             or
           - We hit max_extend expansions from that side.
      5) Return the final best alignment. If we never get good boundaries,
         we keep the best RMSD attempt anyway.

    Args:
        segment_keys (list): The initial missing residues, e.g. [(chain, i, ""), ..., (chain, j, "")]
        filled_res_dict (dict): Residue->atoms for the filled structure
        partial_res_dict (dict): Residue->atoms for the partial structure
        max_extend (int): Max expansions on each side to try.

    Returns:
        (adjusted_atoms, local_rmsd, anchor_atoms)
    """
    import math
    import numpy as np
    from copy import deepcopy

    if not segment_keys:
        raise ValueError("segment_keys is empty.")

    chain = segment_keys[0][0]
    first_missing = segment_keys[0][1]
    last_missing = segment_keys[-1][1]

    print(f"\n=== local_stitch_segment_extend for {chain} {first_missing}-{last_missing} ===")

    # ---------------------------------------------------------------------
    # Step 0: Identify how many anchors we have
    # ---------------------------------------------------------------------
    left_anchor_key = (chain, first_missing - 1, "")
    right_anchor_key = (chain, last_missing + 1, "")
    has_left_anchor = (left_anchor_key in partial_res_dict)
    has_right_anchor = (right_anchor_key in partial_res_dict)
    anchor_count = int(has_left_anchor) + int(has_right_anchor)

    # For returning as anchor_atoms (even if we shift them away)
    left_anchor_atoms = partial_res_dict.get(left_anchor_key, None)
    right_anchor_atoms = partial_res_dict.get(right_anchor_key, None)
    anchor_atoms = (left_anchor_atoms, right_anchor_atoms)

    # If single-res + 2 anchors => use _bridge_single_residue
    if len(segment_keys) == 1 and anchor_count == 2:
        print("  => Single-res bridging with both anchors => _bridge_single_residue")
        single_key = segment_keys[0]
        filled_bb = get_backbone_atoms(filled_res_dict[single_key])
        left_bb = get_backbone_atoms(partial_res_dict[left_anchor_key])
        right_bb = get_backbone_atoms(partial_res_dict[right_anchor_key])
        if ("C" in left_bb and "N" in right_bb and
            "N" in filled_bb and "C" in filled_bb):
            bridged_atoms = _bridge_single_residue(filled_bb, left_bb["C"], right_bb["N"], 1.33)
            return bridged_atoms, 0.0, anchor_atoms

    # If no anchors => copy raw
    if anchor_count == 0:
        print("  => no anchors => copy raw segment atoms.")
        seg_atoms = []
        for key in segment_keys:
            bb_fill = get_backbone_atoms(filled_res_dict[key])
            seg_atoms.extend(list(bb_fill.values()))
        return [a.copy() for a in seg_atoms], 0.0, anchor_atoms

    # If single anchor => do single-boundary approach
    if anchor_count == 1:
        print("  => single-anchor => do minimal alignment + _fix_boundary_bonds")
        # Gather segment backbone
        seg_atoms = []
        for key in segment_keys:
            bb_fill = get_backbone_atoms(filled_res_dict[key])
            seg_atoms.extend(list(bb_fill.values()))
        adjusted = [a.copy() for a in seg_atoms]
        _fix_boundary_bonds(chain, segment_keys, partial_res_dict, adjusted)
        return adjusted, 0.0, anchor_atoms

    # Now, the interesting case: anchor_count == 2 => internal segment
    # We'll do expansions if boundaries are not in range.

    # We'll track how many expansions we've done on each side
    left_expanded = 0
    right_expanded = 0

    def gather_segment_atoms(skeys):
        """Collect the backbone atoms from filled_res_dict for the specified keys."""
        out = []
        for k in skeys:
            bb = get_backbone_atoms(filled_res_dict[k])
            out.extend(list(bb.values()))
        return out

    def build_segment_keys(start_res, end_res):
        """Return a list of keys [ (chain, r, "") for r in range(start_res, end_res+1 ) ]."""
        return [(chain, r, "") for r in range(start_res, end_res+1)]

    def do_kabsch(seg_keys, left_key, right_key):
        """
        Attempt a Kabsch alignment with the two anchor keys (partial->filled).
        If fewer than 2 anchor pairs from those 2 keys => can't do stable rotation.
        Return (adjusted_atoms, local_rmsd, left_dist, right_dist).
        """
        # Build anchor pairs
        anchor_pairs = []
        for akey in [left_key, right_key]:
            if akey in partial_res_dict and akey in filled_res_dict:
                p_bb = get_backbone_atoms(partial_res_dict[akey])
                f_bb = get_backbone_atoms(filled_res_dict[akey])
                for nm in ["N", "CA", "C", "O"]:
                    if nm in p_bb and nm in f_bb:
                        anchor_pairs.append((p_bb[nm], f_bb[nm]))

        seg_atoms_filled = gather_segment_atoms(seg_keys)
        if len(anchor_pairs) < 2:
            # Just copy
            return [a.copy() for a in seg_atoms_filled], 999.9, None, None

        import numpy as np
        partial_coords = []
        filled_coords = []
        for (p_atom, f_atom) in anchor_pairs:
            partial_coords.append(np.array([p_atom["x"], p_atom["y"], p_atom["z"]], dtype=float))
            filled_coords.append(np.array([f_atom["x"], f_atom["y"], f_atom["z"]], dtype=float))
        P = np.array(filled_coords, dtype=float)  # "from"
        Q = np.array(partial_coords, dtype=float) # "to"

        R_local, t_local = kabsch(P, Q)

        # Apply
        adjusted = []
        for atm in seg_atoms_filled:
            new_a = atm.copy()
            c = np.array([atm["x"], atm["y"], atm["z"]])
            new_c = R_local.dot(c) + t_local
            new_a["x"], new_a["y"], new_a["z"] = new_c.tolist()
            adjusted.append(new_a)

        # RMSD
        P_t = (R_local.dot(P.T)).T + t_local
        local_rmsd = float(math.sqrt(np.mean(np.sum((P_t - Q)**2, axis=1))))

        # Check boundary distances
        seg_res = group_by_residue(adjusted)
        left_dist = None
        right_dist = None

        if (chain, seg_keys[0][1], "") in seg_res and left_key in partial_res_dict:
            okL, distL = check_boundary_amide(
                partial_res_dict[left_key],
                seg_res[(chain, seg_keys[0][1], "")]
            )
            left_dist = distL
        if (chain, seg_keys[-1][1], "") in seg_res and right_key in partial_res_dict:
            okR, distR = check_boundary_amide(
                seg_res[(chain, seg_keys[-1][1], "")],
                partial_res_dict[right_key]
            )
            right_dist = distR

        return adjusted, local_rmsd, left_dist, right_dist

    current_first = first_missing
    current_last = last_missing

    best_atoms = []
    best_rmsd = 999.9
    done = False
    iteration = 0

    while iteration < 2 * max_extend + 10:
        iteration += 1
        seg_keys = build_segment_keys(current_first, current_last)
        left_key = (chain, current_first - 1, "")
        right_key = (chain, current_last + 1, "")
        print(f"\nIteration {iteration}: missing={current_first}-{current_last}, anchors={left_key}, {right_key}")

        adj_atoms, local_rmsd, left_dist, right_dist = do_kabsch(seg_keys, left_key, right_key)
        print(f"  => RMSD={local_rmsd:.3f}, left_dist={left_dist}, right_dist={right_dist}")

        if local_rmsd < best_rmsd:
            best_rmsd = local_rmsd
            best_atoms = adj_atoms

        # Evaluate boundary quality
        # If dist is None => maybe the anchor doesn't exist or we had <2 anchor pairs, etc.
        left_good = True
        right_good = True
        if left_dist is not None:
            left_good = (1.2 <= left_dist <= 1.5)
        if right_dist is not None:
            right_good = (1.2 <= right_dist <= 1.5)

        print(f"  => left_good={left_good}, right_good={right_good}")

        if left_dist is None and right_dist is None:
            # We have no valid anchor pairs => break
            print("  => no valid anchors => giving up")
            break

        # If both boundaries are good => done
        if left_good and right_good:
            print("  => Both boundaries good => done")
            best_atoms = adj_atoms
            best_rmsd = local_rmsd
            done = True
            break

        # If left is bad, incorporate the left anchor into missing
        # But only if we haven't extended too far
        if not left_good and left_expanded < max_extend and (left_key in partial_res_dict):
            print(f"  => Incorporating anchor {left_key} into the missing segment on the left")
            current_first = current_first - 1
            left_expanded += 1

        # If right is bad, incorporate the right anchor into missing
        if not right_good and right_expanded < max_extend and (right_key in partial_res_dict):
            print(f"  => Incorporating anchor {right_key} into the missing segment on the right")
            current_last = current_last + 1
            right_expanded += 1

        # If we cannot do anything else => break
        if left_good or left_expanded >= max_extend or (left_key not in partial_res_dict):
            # no more expansions on left side
            pass
        if right_good or right_expanded >= max_extend or (right_key not in partial_res_dict):
            # no more expansions on right side
            pass

        # If we cannot expand either side but boundary not fixed => break
        can_expand_left = (not left_good and left_expanded < max_extend and (left_key in partial_res_dict))
        can_expand_right = (not right_good and right_expanded < max_extend and (right_key in partial_res_dict))
        if not (can_expand_left or can_expand_right):
            print("  => cannot expand further, stopping.")
            break

    if not done:
        print(f"  => Final: boundaries not perfect, best_rmsd={best_rmsd:.3f}, missing={current_first}-{current_last}")

    return best_atoms, best_rmsd, anchor_atoms



def old_local_stitch_segment(segment_keys, filled_res_dict, partial_res_dict, max_range=10):
    """
    Stitch a missing segment from the filled structure into the partial structure
    by iteratively trying anchor residues in expanding ranges from the gap boundaries.

    Workflow:
      1) If exactly one residue is missing and both anchors exist, use _bridge_single_residue.
      2) Otherwise, for range_size in [1..max_range]:
         - Gather partial residues from [first_missing - range_size ... first_missing-1]
           and [last_missing+1 ... last_missing+range_size] that also exist in filled_res_dict.
         - Build anchor pairs for all (N, CA, C, O) found.
         - If <2 anchor pairs, skip this range.
         - Kabsch align the segment to these anchor residues.
         - Check boundary amide distance(s):
            * If there's exactly one anchor, we rely on _fix_boundary_bonds for final
              single-boundary fix.
            * If there are two anchors, we see if the boundary distances are in [1.2..1.5].
         - If acceptable, stop searching further range_sizes.
      3) If no acceptable alignment is found by max_range, use the best attempt.

    Args:
        segment_keys (list): Keys for the missing segment (chain, resSeq, iCode).
        filled_res_dict (dict): Mapping from residue keys -> atoms in the filled structure.
        partial_res_dict (dict): Mapping from residue keys -> atoms in the partial structure.
        max_range (int): Max distance from the segment boundaries to look for anchor residues.

    Returns:
        (adjusted_atoms, local_rmsd, anchor_atoms):
            adjusted_atoms: newly stitched atoms for the segment
            local_rmsd: RMSD of the final alignment to the anchor residues
            anchor_atoms: (left_anchor_atoms, right_anchor_atoms) if they exist
    """
    import math
    from copy import deepcopy
    import numpy as np

    if not segment_keys:
        raise ValueError("No residue keys provided for stitching (segment_keys is empty).")

    chain = segment_keys[0][0]
    first_missing = segment_keys[0][1]
    last_missing = segment_keys[-1][1]

    # Identify direct neighbor anchors (for returning)
    left_anchor_key = (chain, first_missing - 1, "")
    right_anchor_key = (chain, last_missing + 1, "")
    left_anchor_atoms = partial_res_dict.get(left_anchor_key, None)
    right_anchor_atoms = partial_res_dict.get(right_anchor_key, None)
    anchor_atoms = (left_anchor_atoms, right_anchor_atoms)

    # Gather backbone atoms from the filled structure for the missing segment
    seg_atoms = []
    for key in segment_keys:
        bb_filled = get_backbone_atoms(filled_res_dict[key])
        seg_atoms.extend(list(bb_filled.values()))

    print(f"\n=== local_stitch_segment for segment {first_missing}-{last_missing} (chain={chain}) ===")

    # --------------------------------------------------------------
    # 1) Single-res bridging if we have both anchors on a 1-res gap
    # --------------------------------------------------------------
    if len(segment_keys) == 1 and left_anchor_atoms and right_anchor_atoms:
        single_key = segment_keys[0]
        filled_bb = get_backbone_atoms(filled_res_dict[single_key])
        left_bb = get_backbone_atoms(left_anchor_atoms)
        right_bb = get_backbone_atoms(right_anchor_atoms)
        if ("C" in left_bb and "N" in right_bb and
            "N" in filled_bb and "C" in filled_bb):
            print("  Single-res gap + both anchors => _bridge_single_residue()")
            bridged_atoms = _bridge_single_residue(filled_bb, left_bb["C"], right_bb["N"], desired_bond=1.33)
            print("  => single-res bridging complete, forced RMSD=0.0")
            return bridged_atoms, 0.0, anchor_atoms

    # Prepare to store the best alignment in case no perfect boundary match is found
    best_adjusted_atoms = [atom.copy() for atom in seg_atoms]
    best_local_rmsd = float("inf")
    found_good_boundary = False

    # Count how many anchors we have for boundary fix decisions
    has_left_anchor = (left_anchor_key in partial_res_dict)
    has_right_anchor = (right_anchor_key in partial_res_dict)
    anchor_count = sum([has_left_anchor, has_right_anchor])

    # We'll define a helper to check boundary lengths
    def check_boundaries_ok(adjusted_list):
        """
        Return True if both anchors exist and both boundary amides are 1.2..1.5 Å.
        If one anchor, or zero anchors, we don't do a strict check here
        (we rely on _fix_boundary_bonds for single-boundary).
        """
        if anchor_count != 2:
            return True  # skip boundary check if we don't have two anchors
        seg_res = group_by_residue(adjusted_list)
        # boundary left
        left_ok, left_dist = True, None
        if (chain, first_missing, "") in seg_res and has_left_anchor:
            left_ok, left_dist = check_boundary_amide(
                partial_res_dict[left_anchor_key],
                seg_res[(chain, first_missing, "")]
            )
        # boundary right
        right_ok, right_dist = True, None
        if (chain, last_missing, "") in seg_res and has_right_anchor:
            right_ok, right_dist = check_boundary_amide(
                seg_res[(chain, last_missing, "")],
                partial_res_dict[right_anchor_key]
            )

        both_ok = (left_ok and right_ok)
        print(f"    -> Boundaries: left={left_dist}, right={right_dist}, both_ok={both_ok}")
        return both_ok

    # ---------------------------------------------------------------
    # 2) Iteratively expand the anchor search range
    # ---------------------------------------------------------------
    for range_size in range(1, max_range + 1):
        print(f"\n  [range_size={range_size}] Trying anchor residues within +/-{range_size} of gap edges")

        # Collect partial residue keys on the left
        left_candidates = []
        start_left = first_missing - range_size
        end_left = first_missing - 1
        for r in range(start_left, end_left + 1):
            ck = (chain, r, "")
            if ck in partial_res_dict and ck in filled_res_dict:
                left_candidates.append(ck)

        # Collect partial residue keys on the right
        right_candidates = []
        start_right = last_missing + 1
        end_right = last_missing + range_size
        for r in range(start_right, end_right + 1):
            ck = (chain, r, "")
            if ck in partial_res_dict and ck in filled_res_dict:
                right_candidates.append(ck)

        # Combine, build anchor pairs
        all_candidates = list(set(left_candidates + right_candidates))
        all_candidates.sort(key=lambda x: x[1])

        anchor_pairs = []
        for ck in all_candidates:
            p_bb = get_backbone_atoms(partial_res_dict[ck])
            f_bb = get_backbone_atoms(filled_res_dict[ck])
            for nm in ["N", "CA", "C", "O"]:
                if nm in p_bb and nm in f_bb:
                    anchor_pairs.append((p_bb[nm], f_bb[nm]))

        print(f"    -> left_candidates={left_candidates}")
        print(f"    -> right_candidates={right_candidates}")
        print(f"    -> total anchor pairs = {len(anchor_pairs)}")

        # If not enough anchor pairs for stable Kabsch
        if len(anchor_pairs) < 2:
            continue

        # Kabsch alignment
        partial_coords = []
        filled_coords = []
        for (p_atom, f_atom) in anchor_pairs:
            partial_coords.append(np.array([p_atom["x"], p_atom["y"], p_atom["z"]], dtype=float))
            filled_coords.append(np.array([f_atom["x"], f_atom["y"], f_atom["z"]], dtype=float))

        P = np.array(filled_coords, dtype=float)  # "from"
        Q = np.array(partial_coords, dtype=float) # "to"

        R_local, t_local = kabsch(P, Q)

        # Apply transform to seg_atoms
        trial_atoms = []
        for atm in seg_atoms:
            new_atm = atm.copy()
            c = np.array([atm["x"], atm["y"], atm["z"]])
            new_c = R_local.dot(c) + t_local
            new_atm["x"], new_atm["y"], new_atm["z"] = new_c.tolist()
            trial_atoms.append(new_atm)

        # Compute local RMSD
        P_t = (R_local.dot(P.T)).T + t_local
        local_rmsd = float(math.sqrt(np.mean(np.sum((P_t - Q)**2, axis=1))))
        print(f"    -> Kabsch done: local_rmsd={local_rmsd:.3f}")

        # Check boundary distances if we have 2 anchors
        boundaries_ok = check_boundaries_ok(trial_atoms)

        # Update "best" if this RMSD is better or if boundaries are good
        improved = False
        if local_rmsd < best_local_rmsd:
            improved = True
        if improved:
            best_adjusted_atoms = trial_atoms
            best_local_rmsd = local_rmsd

        if boundaries_ok:
            print("    => Boundaries are acceptable, stopping range expansion.")
            found_good_boundary = True
            best_adjusted_atoms = trial_atoms
            best_local_rmsd = local_rmsd
            break

    # End of for range_size

    if not found_good_boundary:
        print(f"  Warning: no anchor range from 1..{max_range} yielded fully acceptable boundaries. "
              "Using best attempt with local_rmsd={:.3f}".format(best_local_rmsd))

    # --------------------------------------------------------------------
    # 3) If there's exactly one anchor boundary, do single-boundary fix
    # --------------------------------------------------------------------
    if anchor_count == 1:
        print("  => single-anchor scenario => calling _fix_boundary_bonds()")
        _fix_boundary_bonds(chain, segment_keys, partial_res_dict, best_adjusted_atoms)
    else:
        print(f"  => anchor_count={anchor_count}, skipping boundary fix")

    return best_adjusted_atoms, best_local_rmsd, anchor_atoms

def transfer_missing_segments(partial_atoms, filled_atoms):
    """
    Transfer missing segments from a filled structure into a partial structure,
    using local_stitch_segment_extend to handle expansions. After local_stitch_segment_extend
    returns the final missing range (which may include some residues that originally existed
    in partial_atoms), we remove the old partial copies of those residues before adding the
    new ones, preventing duplicates.

    Args:
        partial_atoms (list): The partial (experimental) PDB atoms.
        filled_atoms (list): The complete (filled) structure atoms.

    Returns:
        (combined_atoms, transferred_events):
            combined_atoms: Updated list of atom dictionaries
            transferred_events: a list of dicts, each with:
                - 'chain'
                - 'start' (final expanded start residue index)
                - 'end'   (final expanded end residue index)
                - 'rmsd'
                - 'left_extended'
                - 'right_extended'
    """
    # Make an editable copy of the partial atoms
    combined_atoms = list(partial_atoms)
    transferred_events = []

    # Prepare residue dictionaries
    partial_res = group_by_residue(partial_atoms)
    filled_res = group_by_residue(filled_atoms)
    partial_keys = set(partial_res.keys())
    filled_keys = list(filled_res.keys())

    # Identify missing segments
    missing_segments = group_missing_segments(filled_keys, partial_keys)

    for chain, segments in missing_segments.items():
        for segment in segments:
            # Each segment is a list of keys for the missing residues
            partial_res = group_by_residue(combined_atoms)

            # local_stitch_segment_extend now returns up to 7 items:
            # (adjusted_atoms, local_rmsd, anchor_atoms,
            #  final_first, final_last, left_expanded, right_expanded)
            (adjusted_atoms,
             local_rmsd,
             anchor_atoms,
             final_first,
             final_last,
             left_extended,
             right_extended) = local_stitch_segment_extend(
                 segment, filled_res, partial_res, max_extend=10
            )

            # Build the final keys for the newly determined missing range
            final_segment_keys = [
                (chain, r, "") for r in range(final_first, final_last + 1)
            ]

            # -------------------------------------------------
            # Remove old partial atoms for the final_segment_keys
            # so we don't end up with duplicates.
            # -------------------------------------------------
            cleaned_atoms = []
            for atm in combined_atoms:
                # Check if this atom belongs to the newly expanded range
                if (
                    atm["chain_id"] == chain and
                    final_first <= atm["res_seq"] <= final_last
                ):
                    # Skip these: we want to replace them with the new stitched version
                    continue
                cleaned_atoms.append(atm)

            combined_atoms = cleaned_atoms

            # Record event
            event = {
                "chain": chain,
                "start": final_first,
                "end": final_last,
                "rmsd": local_rmsd,
                "left_extended": left_extended,
                "right_extended": right_extended
            }
            transferred_events.append(event)

            # Insert newly stitched atoms for final_segment_keys
            seg_dict = group_by_residue(adjusted_atoms)
            for key in final_segment_keys:
                if key in seg_dict:
                    bb = get_backbone_atoms(seg_dict[key])
                    for name in ["N", "CA", "C", "O"]:
                        if name in bb:
                            new_atom = bb[name].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)

            # Update partial_res after appending
            partial_res = group_by_residue(combined_atoms)

    # Print summary
    print(f"\nInserted {len(transferred_events)} missing segments:")
    for ev in transferred_events:
        ch = ev["chain"]
        st = ev["start"]
        en = ev["end"]
        r = ev["rmsd"]
        le = ev["left_extended"]
        re = ev["right_extended"]
        print(f"  Chain {ch} residues {st}–{en}, local RMSD = {r:.3f} Å")
        if (le > 0) or (re > 0):
            print(f"    => extended left={le}, right={re}")

    return combined_atoms, transferred_events


def most_recent_transfer_missing_segments(partial_atoms, filled_atoms):
    """
    Transfer missing segments from a filled structure into a partial structure,
    using local_stitch_segment_extend to expand or shrink the missing range
    until boundary bonds are acceptable.
    """
    combined_atoms = list(partial_atoms)
    transferred_events = []

    partial_res = group_by_residue(partial_atoms)
    filled_res = group_by_residue(filled_atoms)
    partial_keys = set(partial_res.keys())
    filled_keys = list(filled_res.keys())

    missing_segments = group_missing_segments(filled_keys, partial_keys)

    for chain, segments in missing_segments.items():
        for segment in segments:
            partial_res = group_by_residue(combined_atoms)

            # local_stitch_segment_extend now returns (atoms, rmsd, anchor_atoms,
            # final_first, final_last, left_expanded, right_expanded)
            (adjusted_atoms,
             local_rmsd,
             anchor_atoms,
             final_first,
             final_last,
             left_extended,
             right_extended) = local_stitch_segment_extend(
                 segment, filled_res, partial_res, max_extend=20
            )

            final_segment_keys = [
                (chain, r, "") for r in range(final_first, final_last + 1)
            ]

            # Build an event record
            event = {
                "chain": chain,
                "start": final_first,
                "end": final_last,
                "rmsd": local_rmsd,
                "left_extended": left_extended,
                "right_extended": right_extended
            }
            transferred_events.append(event)

            # Insert backbone atoms for that final range
            seg_dict = group_by_residue(adjusted_atoms)
            for key in final_segment_keys:
                if key in seg_dict:
                    bb = get_backbone_atoms(seg_dict[key])
                    for name in ["N", "CA", "C", "O"]:
                        if name in bb:
                            new_atom = bb[name].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)

            partial_res = group_by_residue(combined_atoms)

    # Now produce the final report with indentation for expansions
    # e.g. "Chain A residues 389–401, local RMSD=0.218" plus lines for expansions
    # We'll just do a print at the end, or you might store it in a string.

    print(f"Inserted {len(transferred_events)} missing segments:")
    for ev in transferred_events:
        ch = ev["chain"]
        st = ev["start"]
        en = ev["end"]
        r = ev["rmsd"]
        le = ev["left_extended"]
        re = ev["right_extended"]
        print(f"  Chain {ch} residues {st}–{en}, local RMSD = {r:.3f} Å")
        if (le > 0) or (re > 0):
            print(f"    => extended left={le}, right={re}")

    return combined_atoms, transferred_events


def oldtransfer_missing_segments(partial_atoms, filled_atoms):
    """
    Transfer missing segments from a filled structure into a partial structure.

    Args:
        partial_atoms (list): partial (experimental) atoms
        filled_atoms (list): filled (complete) structure atoms

    Returns:
        (combined_atoms, transferred_events):
            - combined_atoms: updated list of atoms
            - transferred_events: list of dicts with chain, start, end, rmsd
    """
    combined_atoms = list(partial_atoms)
    transferred_events = []

    # Build residue dictionaries
    partial_res = group_by_residue(partial_atoms)
    filled_res = group_by_residue(filled_atoms)
    partial_keys = set(partial_res.keys())
    filled_keys = list(filled_res.keys())

    # Identify missing segments
    missing_segments = group_missing_segments(filled_keys, partial_keys)

    for chain, segments in missing_segments.items():
        for segment in segments:
            # Recompute partial_res from combined_atoms each time
            partial_res = group_by_residue(combined_atoms)

            # Stitch the segment with the new anchor-based approach
            adjusted_atoms, local_rmsd, _ = local_stitch_segment(
                segment, filled_res, partial_res
            )

            # Record the event
            event = {
                "chain": chain,
                "start": segment[0][1],
                "end": segment[-1][1],
                "rmsd": local_rmsd
            }
            transferred_events.append(event)

            # Add the newly placed backbone atoms
            seg_dict = group_by_residue(adjusted_atoms)
            for key in segment:
                if key in seg_dict:
                    bb = get_backbone_atoms(seg_dict[key])
                    for name in ["N", "CA", "C", "O"]:
                        if name in bb:
                            new_atom = bb[name].copy()
                            new_atom["record_name"] = "ATOM"
                            combined_atoms.append(new_atom)

            # Update partial_res (so next segments see the newly added atoms)
            partial_res = group_by_residue(combined_atoms)

    return combined_atoms, transferred_events


def extract_partial_sequence(res_dict, chain):
    """
    Extract the partial amino acid sequence for a given chain.

    Args:
        res_dict (dict): Mapping of residue keys to lists of atom dicts.
        chain (str): Chain identifier.

    Returns:
        tuple: (sequence (str), sorted_keys (list)) where sequence is the
            one-letter code sequence extracted using three_to_one mapping.
    """
    keys = [k for k in res_dict if k[0] == chain]
    keys = sorted(keys, key=lambda k: (k[1], insertion_order(k[2])))
    seq = ""
    for key in keys:
        res_atoms = res_dict[key]
        res_name = res_atoms[0]["res_name"].upper()
        letter = THREE_TO_ONE.get(res_name, "X")
        seq += letter
    return seq, keys


def find_subsequence_indices(full_seq, subseq):
    """
    Find indices in full_seq that correspond to subseq using dynamic
    programming (Longest Common Subsequence approach).

    Args:
        full_seq (str): The full sequence.
        subseq (str): The subsequence to locate.

    Returns:
        list or None: List of indices (0-indexed) in full_seq that match
            subseq, or None if subseq cannot be fully matched.
    """
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

