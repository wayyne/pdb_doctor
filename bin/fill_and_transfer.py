#!/usr/bin/env python3
"""
fill_and_transfer.py

A script to fill missing residues in a PDB structure and integrate them into
a partial structure using either a local or forge model.
"""

import argparse
import sys
import warnings
import pdbdr as dr
import os

def main():
    """
    The main function to process PDB structures:
    1. Fills missing residues using a selected model.
    2. Performs global alignment and local stitching.
    3. Transfers missing segments into the partial structure.
    4. Validates and writes the final chimeric structure.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Fill missing residues in a PDB structure and transfer them into a "
            "partial structure using either a local or forge model."
        )
    )
    parser.add_argument("partial_pdb", help="Path to a partial PDB file.")
    parser.add_argument("pdb_id", help="PDB ID to fetch the full structure from RCSB.")
    parser.add_argument("output", help="Output PDB file for the final chimeric structure.")
    parser.add_argument(
        "--mode",
        choices=["local", "forge"],
        default="forge",
        help="Model initialization mode: 'local' (default) or 'forge'."
    )
    parser.add_argument(
        "--chain",
        dest="chain_id",
        default=None,
        help="Chain ID from the full PDB to use (others will be stripped). "
             "If omitted, the first ATOM chain in the file will be used."
    )
    args = parser.parse_args()

    print("=" * 80)
    print("Starting the filling, alignment, stitching, and segment transfer process...")

    partial_pdb = args.partial_pdb
    pdb_id = args.pdb_id.strip()
    print(f"Working on {pdb_id}")
    output_pdb = args.output
    mode = args.mode
    chain_id = (args.chain_id or "").strip() if args.chain_id is not None else None
    filled_temp = "filled_temp.pdb"

    # Use the unified fill_structure() function (mode selects local or forge).
    full_seq = dr.fill_structure(mode, pdb_id, filled_temp, chain_id=chain_id)

    try:
        # Load both partial and full (filled) structures.
        partial_atoms = dr.read_pdb(partial_pdb)
        filled_atoms = dr.read_pdb(filled_temp)

        # Ensure the chain IDs are consistent.
        dr.update_chain_id(filled_atoms, partial_atoms[0]["chain_id"])
    except Exception as e:
        print("Error reading PDB files:", e)
        sys.exit(1)

    try:
        # Perform global alignment.
        R, t, common_keys = dr.align_filled_to_partial(partial_atoms, filled_atoms)
    except ValueError as e:
        print("Error in global alignment:", e)
        sys.exit(1)

    # Apply transformation to align the filled structure.
    dr.transform_atoms(filled_atoms, R, t)

    # Compute RMSD to assess alignment quality.
    global_rmsd = dr.compute_alignment_rmsd(partial_atoms, filled_atoms, common_keys)
    print(f"Global alignment RMSD: {global_rmsd:.3f} Å")

    # Transfer missing residues from the filled structure to the partial one.
    combined_atoms, transferred_events = dr.transfer_missing_segments(
        partial_atoms, filled_atoms
    )

    # Check for any amide bond inconsistencies.
    dr.check_amide_bonds(combined_atoms, transferred_events, pdb_id)

    # Save the final chimeric structure.
    dr.write_pdb(output_pdb, combined_atoms)
    print(f"Final chimeric structure written to {output_pdb}")

    os.remove(filled_temp)

if __name__ == "__main__":
    main()

