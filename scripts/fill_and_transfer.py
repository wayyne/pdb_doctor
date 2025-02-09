#!/usr/bin/env python3
"""
fill_and_transfer.py

A script to fill missing residues in a PDB structure and integrate them into
a partial structure using either a local or forge model.
"""

import argparse
import sys
import warnings
import pdbdr as dr  # Custom PDB processing library

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
        default="local",
        help="Model initialization mode: 'local' (default) or 'forge'."
    )
    args = parser.parse_args()

    print("=" * 80)
    print("Starting the filling, alignment, stitching, and segment transfer process...")

    partial_pdb = args.partial_pdb
    pdb_id = args.pdb_id.strip()
    output_pdb = args.output
    mode = args.mode
    filled_temp = "filled_temp.pdb"

    # Use the unified fill_structure() function (mode selects local or forge).
    full_seq = dr.fill_structure(mode, pdb_id, filled_temp)

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

    # Save the transformed (aligned) filled structure.
    fitted_filled_file = "fitted_filled.pdb"
    dr.write_pdb(fitted_filled_file, filled_atoms)
    print(f"Fitted filled structure written to {fitted_filled_file}")

    # Transfer missing residues from the filled structure to the partial one.
    combined_atoms, transferred_events = dr.transfer_missing_segments(
        partial_atoms, filled_atoms
    )

    print(f"Inserted {len(transferred_events)} missing segments:")
    for evt in transferred_events:
        print(
            f"  Chain {evt['chain']} residues {evt['start']}–{evt['end']}, "
            f"local RMSD = {evt['rmsd']:.3f} Å"
        )

    # Check for any amide bond inconsistencies.
    dr.check_amide_bonds(combined_atoms, transferred_events, pdb_id)

    # Check for steric clashes.
    #dr.check_structure_for_clashes(combined_atoms, transferred_events, tolerance=2.25)

    # Save the final chimeric structure.
    dr.write_pdb(output_pdb, combined_atoms)
    print(f"Final chimeric structure written to {output_pdb}")


if __name__ == "__main__":
    # Suppress warnings about missing metadata; we know what we're doing!
    warnings.filterwarnings(
        "ignore",
        message="Entity ID not found in metadata, using None as default"
    )
    main()

