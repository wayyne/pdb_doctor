"""
Module: alignment
Author: [Your Name]
Date: [YYYY-MM-DD]
Description:
    This module provides functions to align protein structures using the 
    Kabsch algorithm. It includes methods to compute centroids, determine
    the optimal rotation and translation between two point sets, transform
    atom coordinates, and compute RMSD between aligned Cα atoms.
"""

import numpy as np
import math

from .pdb_io import group_by_residue, get_backbone_atoms


def compute_centroid(coords):
    """
    Compute the centroid (mean position) of a set of 3D coordinates.

    Args:
        coords (ndarray): An array of shape (N, 3) representing N 3D points.

    Returns:
        ndarray: The centroid of the provided coordinates.
    """
    return np.mean(coords, axis=0)


def kabsch(P, Q):
    """
    Compute the optimal rotation matrix R and translation vector t that
    align two sets of points P and Q using the Kabsch algorithm.

    Args:
        P (ndarray): Array of shape (N, 3) for the first set of points.
        Q (ndarray): Array of shape (N, 3) for the second set of points.

    Returns:
        tuple: A tuple (R, t) where R is a 3x3 rotation matrix and t is a 
            3-element translation vector.
    """
    # Compute centroids of P and Q.
    centroid_P = compute_centroid(P)
    centroid_Q = compute_centroid(Q)

    # Center the points.
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # Compute the covariance matrix.
    C = np.dot(P_centered.T, Q_centered)

    # Perform Singular Value Decomposition.
    V, S, Wt = np.linalg.svd(C)

    # Ensure a right-handed coordinate system.
    d = np.linalg.det(np.dot(Wt.T, V.T))
    if d < 0:
        Wt[-1, :] *= -1

    # Compute rotation matrix and translation vector.
    R = np.dot(Wt.T, V.T)
    t = centroid_Q - np.dot(R, centroid_P)
    return R, t


def transform_atoms(atoms, R, t):
    """
    Apply a rotation and translation to a list of atom dictionaries.

    Each atom dictionary is expected to have numerical keys "x", "y", and "z".

    Args:
        atoms (list): List of atom dictionaries.
        R (ndarray): 3x3 rotation matrix.
        t (ndarray): 3-element translation vector.
    """
    for atom in atoms:
        coord = np.array([atom["x"], atom["y"], atom["z"]])
        new_coord = np.dot(R, coord) + t
        atom["x"], atom["y"], atom["z"] = new_coord.tolist()


def compute_alignment_rmsd(partial_atoms, filled_atoms, common_keys):
    """
    Compute the root-mean-square deviation (RMSD) between the Cα atoms of two
    aligned structures using a set of common residue keys.

    Args:
        partial_atoms (list): List of atom dictionaries from the partial 
            structure.
        filled_atoms (list): List of atom dictionaries from the filled structure.
        common_keys (list): List of residue keys (tuples) common to both 
            structures.

    Returns:
        float: The RMSD value computed over the Cα atoms.
    """
    partial_res = group_by_residue(partial_atoms)
    filled_res = group_by_residue(filled_atoms)
    sq_errors = []
    for key in common_keys:
        pb = get_backbone_atoms(partial_res[key])
        fb = get_backbone_atoms(filled_res[key])
        if "CA" in pb and "CA" in fb:
            diff = (np.array([pb["CA"]["x"], pb["CA"]["y"], pb["CA"]["z"]]) -
                    np.array([fb["CA"]["x"], fb["CA"]["y"], fb["CA"]["z"]]))
            sq_errors.append(np.dot(diff, diff))
    return np.sqrt(np.mean(sq_errors))


def align_filled_to_partial(partial_atoms, filled_atoms):
    """
    Determine the optimal alignment of a filled structure to a partial structure
    based on common Cα atoms using the Kabsch algorithm.

    Args:
        partial_atoms (list): List of atom dictionaries for the partial structure.
        filled_atoms (list): List of atom dictionaries for the filled structure.

    Returns:
        tuple: A tuple (R, t, common_keys) where:
            - R is the 3x3 rotation matrix.
            - t is the 3-element translation vector.
            - common_keys is a list of residue keys used for the alignment.
    """
    partial_res = group_by_residue(partial_atoms)
    filled_res = group_by_residue(filled_atoms)
    common_keys = []
    partial_coords = []
    filled_coords = []

    # Loop over residues in the partial structure and identify common ones.
    for key in partial_res:
        if key in filled_res:
            pb = get_backbone_atoms(partial_res[key])
            fb = get_backbone_atoms(filled_res[key])
            if "CA" in pb and "CA" in fb:
                common_keys.append(key)
                partial_coords.append(np.array(
                    [pb["CA"]["x"], pb["CA"]["y"], pb["CA"]["z"]]
                ))
                filled_coords.append(np.array(
                    [fb["CA"]["x"], fb["CA"]["y"], fb["CA"]["z"]]
                ))
    # Ensure that there are enough common points for alignment.
    if len(partial_coords) < 3:
        raise ValueError("Not enough common CA atoms for global alignment.")

    P = np.array(filled_coords)
    Q = np.array(partial_coords)
    R, t = kabsch(P, Q)
    return R, t, common_keys

