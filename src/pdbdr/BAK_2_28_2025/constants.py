"""
Module: constants
Author: [Your Name]
Date: [YYYY-MM-DD]
Description:
    This module centralizes all constant definitions for the pdbdr package.
    It includes API endpoints, atom type lists and mappings, amino acid code
    maps, backbone atom names, van der Waals radii, and standard bond lengths
    and angles for peptide structures.
"""

import math

# --- Forge API Constants ---
FORGE_URL = "https://forge.evolutionaryscale.ai"
FORGE_MDL = "esm3-large-2024-03"

# --- Atom Types and Ordering ---
ATOM_TYPES = [
    "N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1", "SG",
    "CD", "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1",
    "CE2", "CE3", "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2",
    "OH", "CZ", "CZ2", "CZ3", "NZ", "OXT"
]
ATOM_ORDER = {atom: i for i, atom in enumerate(ATOM_TYPES)}
ATOM_TYPE_NUM = len(ATOM_TYPES)

# --- Amino Acid Code Mappings ---
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")  # Set of valid amino acid characters

ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
    "U": "SEC",  # Selenocysteine
    "O": "PYL",  # Pyrrolysine
    "X": "UNK",  # Unknown residue
    "-": "UNK"   # Explicitly handle gaps
}

# --- Backbone and Structural Constants ---
BACKBONE_NAMES = ["N", "CA", "C", "O"]

# Van der Waals radii (in Å) for common elements (used in clash detection)
VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
}

# Standard bond lengths (in Å) and bond angles (in radians) used for
# building ideal peptide fragments.
BOND_N_CA = 1.46
BOND_CA_C = 1.52
BOND_C_O = 1.23

ANGLE_NCAC = math.radians(111.0)
ANGLE_CCO = math.radians(120.0)

