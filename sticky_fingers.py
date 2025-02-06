#!/usr/bin/env python3
import os
import glob
import math
import sys

##########################################
# Reused function: read_pdb (provided)
##########################################
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

##########################################
# Helper: write_pdb
##########################################
def write_pdb(filename, atoms):
    """Write a list of atoms to a PDB file."""
    with open(filename, "w") as f:
        for atom in atoms:
            f.write(
                f"{atom['record_name']:6}{atom['serial']:5} {atom['atom_name']:<4}{atom['alt_loc']:1}"
                f"{atom['res_name']:3} {atom['chain_id']:1}{atom['res_seq']:4}{atom['i_code']:1}   "
                f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}{atom['occupancy']:6.2f}{atom['temp_factor']:6.2f}          "
                f"{atom['element']:>2}{atom['charge']:>2}\n"
            )

##########################################
# Step 1: Transfer HETATM records
##########################################
def transfer_hetatm_records(output_dir, og_pp_pdbs_dir):
    """
    Transfers HETATM records from the original pickpocket PDBs (in og_pp_pdbs_dir)
    to their counterparts in output_dir.
    """
    for pdb_file in os.listdir(output_dir):
        if pdb_file.endswith("_C.pdb"):
            base_name = pdb_file[:7]  # Extract the first 7 characters
            og_pdb_path = os.path.join(og_pp_pdbs_dir, f"{base_name}.pdb")
            output_pdb_path = os.path.join(output_dir, pdb_file)

            if os.path.exists(og_pdb_path):
                hetatm_records = []
                with open(og_pdb_path, "r") as f:
                    for line in f:
                        if line.startswith("HETATM"):
                            hetatm_records.append(line)
                # Append the HETATM records to the output file
                with open(output_pdb_path, "a") as f:
                    f.writelines(hetatm_records)
                print(f"Transferred HETATM records from {og_pdb_path} to {output_pdb_path}")

##########################################
# Step 2: Compute Contacts and Relabel PDB
##########################################
def compute_contacts(atoms, contact_threshold=5.0):
    """
    Identify protein residues in contact with any ligand heavy atom.
    A contact is defined as any heavy atom (non-hydrogen) in a ligand (HETATM)
    being within contact_threshold Å of any heavy atom in a residue (ATOM record).
    
    Returns a set of residue identifiers (tuple of chain_id, res_seq, i_code) that are in contact.
    """
    def distance(atom1, atom2):
        return math.sqrt(
            (atom1["x"] - atom2["x"]) ** 2 +
            (atom1["y"] - atom2["y"]) ** 2 +
            (atom1["z"] - atom2["z"]) ** 2
        )

    # Consider only heavy atoms (skip hydrogens) -- JUST KIDDING FOR NOW
    #ligand_atoms = [atom for atom in atoms if atom["record_name"] == "HETATM" and atom["element"].upper() != "H"]
    #protein_atoms = [atom for atom in atoms if atom["record_name"] == "ATOM" and atom["element"].upper() != "H"]
    ligand_atoms = [atom for atom in atoms if atom["record_name"] == "HETATM"]
    protein_atoms = [atom for atom in atoms if atom["record_name"] == "ATOM"]

    contact_residues = set()
    for latom in ligand_atoms:
        for patom in protein_atoms:
            if distance(latom, patom) < contact_threshold:
                # Use (chain_id, res_seq, i_code) as the unique identifier for a residue.
                contact_residues.add((patom["chain_id"], patom["res_seq"], patom["i_code"]))
    return contact_residues

def relabel_pdb_with_contacts(pdb_path):
    """
    Reads a PDB file, computes contacts, and updates beta factors:
      - Sets temp_factor to 1.0 for ATOM records belonging to residues in contact.
      - Sets temp_factor to 0.0 for ATOM records not in contact.
    The updated file is written back to pdb_path.
    """
    atoms = read_pdb(pdb_path)
    contact_residues = compute_contacts(atoms, contact_threshold=5.0)

    for atom in atoms:
        if atom["record_name"] == "ATOM":
            key = (atom["chain_id"], atom["res_seq"], atom["i_code"])
            if key in contact_residues:
                atom["temp_factor"] = 1.0
            else:
                atom["temp_factor"] = 0.0

    write_pdb(pdb_path, atoms)
    print(f"Updated {pdb_path} with contact relabeling.")

##########################################
# Helper: Compute minimum distance for a residue
##########################################
def min_distance_for_residue(atoms, residue_key):
    """
    For a given residue (identified by (chain_id, res_seq, i_code)), compute
    the minimum distance between any heavy atom (non-hydrogen) in that residue
    (from ATOM records) and any heavy atom in the ligand (HETATM records).
    
    Returns the minimum distance as a float, or None if no atoms found.
    """
    residue_atoms = [
        atom for atom in atoms
        if atom["record_name"] == "ATOM" 
        and (atom["chain_id"], atom["res_seq"], atom["i_code"]) == residue_key
        and atom["element"].upper() != "H"
    ]
    ligand_atoms = [
        atom for atom in atoms
        if atom["record_name"] == "HETATM"
        and atom["element"].upper() != "H"
    ]
    if not residue_atoms or not ligand_atoms:
        return None
    min_dist = float('inf')
    for latom in ligand_atoms:
        for patom in residue_atoms:
            dx = latom["x"] - patom["x"]
            dy = latom["y"] - patom["y"]
            dz = latom["z"] - patom["z"]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < min_dist:
                min_dist = dist
    return min_dist

##########################################
# Helper: Count statistics from a list of atoms
##########################################
def count_stats(atoms):
    """
    Returns a tuple with:
      - heavy_atoms: count of atoms whose element is not "H"
      - residue_count: count of unique amino acid residues (from ATOM records) 
                       identified by (chain_id, res_seq, i_code)
      - positive_labels: count of residues (from ATOM records) with temp_factor == 1.0
      - negative_labels: count of residues (from ATOM records) with temp_factor == 0.0
    """
    heavy_atoms = sum(1 for atom in atoms if atom["element"].upper() != "H")
    residue_dict = {}
    for atom in atoms:
        if atom["record_name"] == "ATOM":
            key = (atom["chain_id"], atom["res_seq"], atom["i_code"])
            if key not in residue_dict:
                residue_dict[key] = atom["temp_factor"]
    residue_count = len(residue_dict)
    positive_labels = sum(1 for v in residue_dict.values() if v == 1.0)
    negative_labels = sum(1 for v in residue_dict.values() if v == 0.0)
    return heavy_atoms, residue_count, positive_labels, negative_labels

##########################################
# Helper: Extract positive residue set from ATOM records
##########################################
def get_positive_residues(atoms):
    """
    Returns a set of residue identifiers (chain_id, res_seq, i_code) from ATOM records
    that have temp_factor == 1.0.
    """
    pos_set = {}
    for atom in atoms:
        if atom["record_name"] == "ATOM":
            key = (atom["chain_id"], atom["res_seq"], atom["i_code"])
            # If we haven't seen the residue, record its temp_factor.
            if key not in pos_set:
                pos_set[key] = atom["temp_factor"]
    return {key for key, val in pos_set.items() if val == 1.0}

##########################################
# Main: Process each output PDB and log details
##########################################
def main():
    output_dir = sys.argv[1]
    if len(sys.argv) < 2:
        print("Usage: sticky_fingers.py <pdb_dir>")
        sys.exit(1) 

    output_dir = sys.argv[1] 
    print(f"Output directory: {output_dir}")

    basecamp = "/home/wayyne/wrk/pdb_doctor"
    og_pp_pdbs_dir = "/home/wayyne/wrk/biolip2_for_pretraining/pickpocket_2025/pp_pdbs/single-chain"

    # Open log files
    log_filename = "sticky_fingers_log.txt"
    diff_filename = "sticky_fingers_diff.txt"
    log_file = open(log_filename, "w")
    diff_file = open(diff_filename, "w")

    # Write header for diff.txt
    diff_file.write("Filename\tHeavyAtomsDiff\tResiduesDiff\tPositivesDiff\tNegativesDiff\n")

    # Step 1: Transfer HETATM records
    transfer_hetatm_records(output_dir, og_pp_pdbs_dir)

    # Process each output PDB file
    for pdb_file in os.listdir(output_dir):
        if pdb_file.endswith("_C.pdb"):
            output_pdb_path = os.path.join(output_dir, pdb_file)
            base_name = pdb_file[:4]
            reference_pdb_path = glob.glob(os.path.join(og_pp_pdbs_dir, f"{base_name}*.pdb"))[0]
            #reference_pdb_path = os.path.join(og_pp_pdbs_dir, f"{base_name}*.pdb")

            log_file.write("="*50 + "\n")
            log_file.write(f"Processing file: {pdb_file}\n")
            log_file.write(f"Reference PDB: {reference_pdb_path}\n")

            # Read reference PDB stats (if exists)
            if os.path.exists(reference_pdb_path):
                ref_atoms = read_pdb(reference_pdb_path)
                ref_stats = count_stats(ref_atoms)
                # Get positive residues from the reference file (based on stored beta-factors)
                ref_positive_set = get_positive_residues(ref_atoms)
            else:
                log_file.write("Reference PDB not found. Skipping stats for reference.\n")
                ref_stats = (0, 0, 0, 0)
                ref_positive_set = set()

            # Relabel output PDB (this updates the file on disk)
            relabel_pdb_with_contacts(output_pdb_path)

            # Re-read updated output PDB
            out_atoms = read_pdb(output_pdb_path)
            out_stats = count_stats(out_atoms)
            # Get positive residues from the new file (based on stored beta-factors)
            new_positive_set = get_positive_residues(out_atoms)

            # Compute differences (final - reference)
            diff_heavy = out_stats[0] - ref_stats[0]
            diff_residues = out_stats[1] - ref_stats[1]
            diff_positives = out_stats[2] - ref_stats[2]
            diff_negatives = out_stats[3] - ref_stats[3]

            # Write details to the log file
            log_file.write(f"Heavy Atoms: Reference = {ref_stats[0]}, Final = {out_stats[0]}, Diff = {diff_heavy}\n")
            log_file.write(f"Residues:    Reference = {ref_stats[1]}, Final = {out_stats[1]}, Diff = {diff_residues}\n")
            log_file.write(f"Positives:   Reference = {ref_stats[2]}, Final = {out_stats[2]}, Diff = {diff_positives}\n")
            log_file.write(f"Negatives:   Reference = {ref_stats[3]}, Final = {out_stats[3]}, Diff = {diff_negatives}\n")

            # Identify residues that were positive in the reference but not in the new file
            lost_positives = ref_positive_set - new_positive_set
            if lost_positives:
                log_file.write("\nResidues positive in reference but not in new PDB:\n")
                for res_key in lost_positives:
                    ref_min = min_distance_for_residue(ref_atoms, res_key)
                    new_min = min_distance_for_residue(out_atoms, res_key)
                    if ref_min is not None and new_min is not None:
                    		log_file.write(f"Residue {res_key}: Ref min distance = {ref_min:.3f} Å, New min distance = {new_min:.3f} Å\n")
                    else:
                    		log_file.write(f"Residue {res_key}: Distance calculation failed (None value present)\n")
 
            else:
                log_file.write("\nNo lost positive residues detected.\n")

            log_file.write("="*50 + "\n\n")

            # Write a summary line to diff.txt using formatted columns
            diff_file.write(f"{base_name}\t{diff_heavy:5d}\t{diff_residues:5d}\t{diff_positives:5d}\t{diff_negatives:5d}\n")

    log_file.close()
    diff_file.close()
    print(f"Processing complete. See {log_filename} and {diff_filename} for details.")

if __name__ == "__main__":
    main()

