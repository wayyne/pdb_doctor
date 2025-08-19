#!/usr/bin/env python3

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
    "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
    "VAL": "V",
}

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

ONE_TO_THREE = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    "U": "SEC",  # Selenocysteine
    "O": "PYL",  # Pyrrolysine
    "X": "UNK",  # Unknown residue
    "-": "UNK"
}

# --- Modified Residues Mapping ---
# --- Modified Residues Mapping ---
MODIFIED_RESIDUE_MAP = {
    # Phosphorylations
    "SEP": "S",  # Phosphoserine
    "TPO": "T",  # Phosphothreonine
    "PTR": "Y",  # Phosphotyrosine

    # Oxidations and Hydroxylations
    "CSO": "C",  # S-hydroxycysteine
    "HYP": "P",  # Hydroxyproline
    "MLY": "K",  # N6-methyllysine
    "M3L": "K",  # N6,N6,N6-trimethyllysine
    "HIC": "H",  # 4-methyl-histidine

    # Seleno and sulfur analogs
    "MSE": "M",  # Selenomethionine
    "SEC": "C",  # Selenocysteine (map to cysteine if ESM can't handle SEC)

    # Modified cysteines
    "CME": "C",  # S,S-(2-hydroxyethyl)cysteine
    "CSX": "C",  # S-oxycysteine
    "CSS": "C",  # Cystine
    "ACY": "C",  # S-acetylcysteine

    # Acetylations and methylations
    "ALY": "K",  # N-acetyllysine
    "MLZ": "K",  # Methyllysine
    "M2L": "K",  # Dimethyllysine
    "M3L": "K",  # Trimethyllysine
    "AIB": "A",  # Alpha-aminoisobutyric acid (closest: alanine)

    # Glycosylations (map to base residue)
#    "NAG": "N",  # N-acetylglucosamine-modified asparagine (approx. map)
    "BMA": "N",  # Beta-mannose-modified Asn

    # Formyl, pyro-Glu, and others
    "FME": "M",  # N-formylmethionine
    "PCA": "E",  # Pyro-glutamic acid from glutamate
    "OCS": "C",  # Cysteine sulfinic acid

    # Miscellaneous
    "ORN": "R",  # Ornithine (nonstandard, closest: arginine)
    "DPR": "P",  # D-proline (assume normal L-proline)
    "DGL": "E",  # D-glutamate (map to E)
    "DGN": "Q",  # D-glutamine (map to Q)
}

# --- Backbone and Structural Constants ---
BACKBONE_NAMES = ["N", "CA", "C", "O"]

VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75,
}

BOND_N_CA = 1.46
BOND_CA_C = 1.52
BOND_C_O = 1.23

ANGLE_NCAC = math.radians(111.0)
ANGLE_CCO = math.radians(120.0)

