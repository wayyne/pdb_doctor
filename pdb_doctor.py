#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from collections import defaultdict
import shutil
import concurrent.futures
import logging
import statistics

# Macros for output suffixes
OUTPUT_A_SUFFIX = "A"
OUTPUT_B_SUFFIX = "B"
OUTPUT_C_SUFFIX = "C"

def main():
    parser = argparse.ArgumentParser(
        description=("Process PDB files: add missing heavy atoms, optimize with "
                     "SCWRL4, and merge coordinates, or restore complete missing "
                     "residues.")
    )
    parser.add_argument('--input', required=True,
                        help="Input directory containing PDB files")
    parser.add_argument('--output', required=True,
                        help="Output directory to save processed PDB files")
    parser.add_argument(
        '--mode',
        choices=['partial', 'complete'],
        default='complete',
        help=("Select restoration mode: 'partial' to fix partial missing "
              "electron density; 'complete' to restore complete missing "
              "residues (default: complete)"))
    parser.add_argument(
        '--restore-mode',
        choices=['termini', 'gaps', 'both'],
        default='both',
        help=("For complete mode, select which missing residues to restore: "
              "termini, gaps, or both (default: both)"))
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.output

    # Ensure SCWRL4 is in PATH
    scwrl4_path = shutil.which('Scwrl4')
    if scwrl4_path is None:
        print("Error: SCWRL4 executable 'Scwrl4' not found in PATH.")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # Set up logging
    log_file = os.path.join(output_dir, "fix_pdbs.log")
    setup_logging(log_file)
    logger = logging.getLogger()

    # List all PDB files in the input directory
    pdb_files = [f for f in os.listdir(input_dir)
                 if f.lower().endswith('.pdb')]
    if not pdb_files:
        logger.info(f"No PDB files found in the input directory: {input_dir}")
        sys.exit(1)

    logger.info(f"Found {len(pdb_files)} PDB files in '{input_dir}'. "
                "Starting processing...\n")

    # Statistics
    total_files = len(pdb_files)
    success_count = 0
    failure_count = 0
    missing_atoms_counts = []

    # Process PDB files in parallel
    with concurrent.futures.ProcessPoolExecutor() as executor:
        if args.mode == 'partial':
            futures = {executor.submit(process_pdb_file, pdb_file, input_dir,
                                         output_dir, scwrl4_path): pdb_file
                       for pdb_file in pdb_files}
        else:  # complete mode
            futures = {executor.submit(process_pdb_file_complete, pdb_file,
                                         input_dir, output_dir, scwrl4_path,
                                         args.restore_mode): pdb_file
                       for pdb_file in pdb_files}
        for future in concurrent.futures.as_completed(futures):
            pdb_file = futures[future]
            try:
                success, log_messages, missing_atom_count = future.result()
                logger.info("\n".join(log_messages))
                if success:
                    success_count += 1
                    missing_atoms_counts.append(missing_atom_count)
                else:
                    failure_count += 1
            except Exception as e:
                failure_count += 1
                logger.error(f"Error processing {pdb_file}: {e}")

    # Calculate statistics
    if missing_atoms_counts:
        average_missing_atoms = statistics.mean(missing_atoms_counts)
        if len(missing_atoms_counts) > 1:
            std_dev_missing_atoms = statistics.stdev(missing_atoms_counts)
        else:
            std_dev_missing_atoms = 0.0
    else:
        average_missing_atoms = 0.0
        std_dev_missing_atoms = 0.0

    # Summary
    summary = (
        "\nProcessing Complete.\n"
        f"  Total files processed: {total_files}\n"
        f"  Successfully processed: {success_count}\n"
        f"  Failed to process: {failure_count}\n"
        f"  Average number of missing atoms per file: "
        f"{average_missing_atoms:.2f}\n"
        f"  Standard deviation of missing atoms per file: "
        f"{std_dev_missing_atoms:.2f}\n"
    )

    print(summary)
    logger.info(summary)
    logger.info(f"Detailed logs are available in '{log_file}'.")
    # Optionally, write statistics to a separate 'stats.dat' file
    stats_file = os.path.join(output_dir, "stats.dat")
    with open(stats_file, 'w') as f:
        f.write(summary)
    logger.info(f"Statistics have been written to '{stats_file}'.")

def process_pdb_file(pdb_file, input_dir, output_dir, scwrl4_path):
    log_messages = []
    input_pdb_path = os.path.join(input_dir, pdb_file)
    base_name = os.path.splitext(pdb_file)[0]
    output_pdb_a = os.path.join(output_dir,
                                f"{base_name}_{OUTPUT_A_SUFFIX}.pdb")
    output_pdb_b = os.path.join(output_dir,
                                f"{base_name}_{OUTPUT_B_SUFFIX}.pdb")
    output_pdb_c = os.path.join(output_dir,
                                f"{base_name}_{OUTPUT_C_SUFFIX}.pdb")

    log_messages.append(f"Processing {pdb_file}...")

    try:
        # Read original input PDB atoms with beta-factors
        original_atoms = read_pdb_atoms(input_pdb_path)
        beta_factors, beta_factors_residue = get_beta_factors(original_atoms)

        # Step 1: Run PDBFixer to add missing heavy atoms
        fixer = PDBFixer(filename=input_pdb_path)
        fixer.findMissingResidues()
        if fixer.missingResidues:
            log_messages.append("  Missing residues detected:")
            for chain_id in fixer.missingResidues:
                res_ids = [str(res_id) for res_id in fixer.missingResidues[chain_id]]
                log_messages.append(f"    Chain {chain_id}: Residues "
                                    f"{', '.join(res_ids)}")
        else:
            log_messages.append("  No missing residues detected.")

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        if fixer.missingAtoms:
            missing_atoms_total = 0
            log_messages.append(f"  Missing atoms detected for "
                                f"{len(fixer.missingAtoms)} residues:")
            for residue, atoms in fixer.missingAtoms.items():
                chain_id = residue.chain.id
                res_seq = residue.id
                res_name = residue.name
                atom_names = [atom.name for atom in atoms]
                missing_atoms_total += len(atom_names)
                log_messages.append(f"    Residue {res_name} {res_seq} "
                                    f"(Chain {chain_id}): Missing atoms "
                                    f"{', '.join(atom_names)}")
        else:
            missing_atoms_total = 0
            log_messages.append("  No missing atoms detected.")

        # Do not add hydrogens
        # fixer.addMissingHydrogens()
        # Add missing heavy atoms only
        fixer.addMissingAtoms()

        # Save PDBFixer output (outputPDB-A)
        with open(output_pdb_a, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

        # Transfer beta-factors to outputPDB-A
        fixed_atoms = read_pdb_atoms(output_pdb_a)
        fixed_atoms = assign_beta_factors(fixed_atoms, beta_factors,
                                          beta_factors_residue)
        write_pdb_atoms(fixed_atoms, output_pdb_a)

        log_messages.append(f"  Missing heavy atoms added and saved to "
                            f"{output_pdb_a}")

        # Step 2: Run SCWRL4 on outputPDB-A to get outputPDB-B
        scwrl_command = ['Scwrl4', '-i', output_pdb_a, '-o', output_pdb_b]
        result = subprocess.run(scwrl_command, capture_output=True,
                                text=True)
        if result.returncode != 0:
            log_messages.append(f"Error running SCWRL4 on {output_pdb_a}:\n"
                                f"{result.stderr}")
            return False, log_messages, missing_atoms_total

        # Transfer beta-factors to outputPDB-B
        scwrl_atoms = read_pdb_atoms(output_pdb_b)
        scwrl_atoms = assign_beta_factors(scwrl_atoms, beta_factors,
                                          beta_factors_residue)
        write_pdb_atoms(scwrl_atoms, output_pdb_b)

        log_messages.append(f"  SCWRL4 optimization completed and saved to "
                            f"{output_pdb_b}")

        # Step 3: Create outputPDB-C by merging coordinates and transferring
        # beta-factors
        create_output_pdb_c(original_atoms, fixed_atoms, scwrl_atoms,
                            output_pdb_c)
        log_messages.append(f"  Merged PDB saved to {output_pdb_c}")

        return True, log_messages, missing_atoms_total

    except Exception as e:
        log_messages.append(f"Error processing {pdb_file}: {e}")
        return False, log_messages, 0

def process_pdb_file_complete(pdb_file, input_dir, output_dir, 
                              scwrl4_path, restore_mode):
    log_messages = []
    input_pdb_path = os.path.join(input_dir, pdb_file)
    base_name = os.path.splitext(pdb_file)[0]
    output_pdb_a = os.path.join(output_dir,
                                f"{base_name}_{OUTPUT_A_SUFFIX}_complete.pdb")
    output_pdb_b = os.path.join(output_dir,
                                f"{base_name}_{OUTPUT_B_SUFFIX}_complete.pdb")
    output_pdb_c = os.path.join(output_dir,
                                f"{base_name}_{OUTPUT_C_SUFFIX}_complete.pdb")
    log_messages.append(f"Processing {pdb_file} for complete missing "
                        "residues...")

    try:
        # Read original input PDB atoms with beta-factors
        original_atoms = read_pdb_atoms(input_pdb_path)
        beta_factors, beta_factors_residue = get_beta_factors(original_atoms)

        # Step 1: Run PDBFixer to restore missing residues and add heavy atoms
        fixer = PDBFixer(filename=input_pdb_path)
        fixer.findMissingResidues()

        # Filter missing residues based on restore_mode
        filtered_missing = {}
        # Iterate over chains in the current structure
        for chain in fixer.topology.chains():
            chain_id = chain.id
            present_res = []
            for residue in chain.residues():
                try:
                    # residue.id is a tuple: (het, resSeq, iCode)
                    rseq = int(residue.id[1])
                    present_res.append(rseq)
                except Exception:
                    continue
            if not present_res:
                continue
            min_res = min(present_res)
            max_res = max(present_res)
            if chain_id in fixer.missingResidues:
                missing_res = fixer.missingResidues[chain_id]
                filtered = []
                for res in missing_res:
                    try:
                        res_int = int(res)
                    except Exception:
                        continue
                    if restore_mode == 'termini':
                        if res_int < min_res or res_int > max_res:
                            filtered.append(res)
                    elif restore_mode == 'gaps':
                        if min_res < res_int < max_res:
                            filtered.append(res)
                    else:  # both
                        filtered.append(res)
                if filtered:
                    filtered_missing[chain_id] = filtered

        fixer.missingResidues = filtered_missing

        if filtered_missing:
            log_messages.append("  Missing residues to be restored:")
            for chain_id in filtered_missing:
                res_ids = [str(r) for r in filtered_missing[chain_id]]
                log_messages.append(f"    Chain {chain_id}: Residues "
                                    f"{', '.join(res_ids)}")
        else:
            log_messages.append("  No missing residues to restore based on "
                                f"restore mode '{restore_mode}'.")

        # Insert missing residues. This assumes PDBFixer has this method.
        fixer.insertMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        # Save PDBFixer output (outputPDB-A)
        with open(output_pdb_a, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

        # Transfer beta-factors to outputPDB-A
        fixed_atoms = read_pdb_atoms(output_pdb_a)
        fixed_atoms = assign_beta_factors(fixed_atoms, beta_factors,
                                          beta_factors_residue)
        write_pdb_atoms(fixed_atoms, output_pdb_a)

        log_messages.append(f"  Missing residues and heavy atoms added and "
                            f"saved to {output_pdb_a}")

        # Step 2: Run SCWRL4 on outputPDB-A to get outputPDB-B
        scwrl_command = ['Scwrl4', '-i', output_pdb_a, '-o', output_pdb_b]
        result = subprocess.run(scwrl_command, capture_output=True,
                                text=True)
        if result.returncode != 0:
            log_messages.append(f"Error running SCWRL4 on {output_pdb_a}:\n"
                                f"{result.stderr}")
            return False, log_messages, 0

        scwrl_atoms = read_pdb_atoms(output_pdb_b)
        scwrl_atoms = assign_beta_factors(scwrl_atoms, beta_factors,
                                          beta_factors_residue)
        write_pdb_atoms(scwrl_atoms, output_pdb_b)
        log_messages.append(f"  SCWRL4 optimization completed and saved to "
                            f"{output_pdb_b}")

        # Step 3: Create outputPDB-C by merging coordinates and transferring
        # beta-factors
        create_output_pdb_c(original_atoms, fixed_atoms, scwrl_atoms,
                            output_pdb_c)
        log_messages.append(f"  Merged PDB saved to {output_pdb_c}")

        # For complete mode, missing_atoms_total is not directly applicable.
        return True, log_messages, 0

    except Exception as e:
        log_messages.append(f"Error processing {pdb_file}: {e}")
        return False, log_messages, 0

def get_beta_factors(atoms):
    beta_factors = {}
    beta_factors_residue = defaultdict(list)
    for atom in atoms:
        atom_key = (atom['chain_id'], atom['res_seq'],
                    atom['i_code'], atom['atom_name'])
        beta_factors[atom_key] = atom['temp_factor']
        residue_key = (atom['chain_id'], atom['res_seq'], atom['i_code'])
        beta_factors_residue[residue_key].append(atom['temp_factor'])
    return beta_factors, beta_factors_residue

def assign_beta_factors(atoms, beta_factors, beta_factors_residue):
    for atom in atoms:
        atom_key = (atom['chain_id'], atom['res_seq'],
                    atom['i_code'], atom['atom_name'])
        residue_key = (atom['chain_id'], atom['res_seq'], atom['i_code'])
        if atom_key in beta_factors:
            beta = beta_factors[atom_key]
        elif beta_factors_residue[residue_key]:
            beta = sum(beta_factors_residue[residue_key]) / \
                   len(beta_factors_residue[residue_key])
        else:
            beta = 0.00
        atom['temp_factor'] = beta
    return atoms

def create_output_pdb_c(original_atoms, fixed_atoms, scwrl_atoms,
                        output_pdb_c):
    # Identify atoms that were missing in the original input PDB
    original_atom_keys = {(atom['chain_id'], atom['res_seq'],
                           atom['i_code'], atom['atom_name'])
                          for atom in original_atoms}

    # Create a mapping from atom key to atom for SCWRL atoms
    scwrl_atom_dict = {(atom['chain_id'], atom['res_seq'],
                        atom['i_code'], atom['atom_name']): atom
                       for atom in scwrl_atoms}

    # Create a mapping of beta-factors from the original PDB
    beta_factors, beta_factors_residue = get_beta_factors(original_atoms)

    # Create outputPDB-C by merging coordinates and transferring beta-factors
    output_atoms = []
    for atom in fixed_atoms:
        atom_key = (atom['chain_id'], atom['res_seq'],
                    atom['i_code'], atom['atom_name'])
        residue_key = (atom['chain_id'], atom['res_seq'], atom['i_code'])

        if atom_key in original_atom_keys:
            # Atom was present in the original input PDB; use coordinates from
            # fixed_atoms and transfer beta-factor from original atom.
            atom['temp_factor'] = beta_factors.get(atom_key, 0.00)
            output_atoms.append(atom)
        else:
            # Atom was missing in the original input PDB; use coordinates from
            # SCWRL4 output.
            if atom_key in scwrl_atom_dict:
                new_atom = scwrl_atom_dict[atom_key]
                if beta_factors_residue[residue_key]:
                    avg_beta = sum(beta_factors_residue[residue_key]) / \
                               len(beta_factors_residue[residue_key])
                else:
                    avg_beta = 0.00
                new_atom['temp_factor'] = avg_beta
                output_atoms.append(new_atom)
            else:
                if beta_factors_residue[residue_key]:
                    avg_beta = sum(beta_factors_residue[residue_key]) / \
                               len(beta_factors_residue[residue_key])
                else:
                    avg_beta = 0.00
                atom['temp_factor'] = avg_beta
                output_atoms.append(atom)

    write_pdb_atoms(output_atoms, output_pdb_c)

def read_pdb_atoms(pdb_file_path):
    atoms = []
    with open(pdb_file_path, 'r') as pdb_file:
        serial_counter = 1  # To handle missing serial numbers
        for line in pdb_file:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    serial = line[6:11].strip()
                    serial = int(serial) if serial else serial_counter
                    atom = {
                        'record_name': line[0:6].strip(),
                        'serial': serial,
                        'atom_name': line[12:16].strip(),
                        'alt_loc': line[16].strip(),
                        'res_name': line[17:20].strip(),
                        'chain_id': line[21].strip(),
                        'res_seq': line[22:26].strip(),
                        'i_code': line[26].strip(),
                        'x': float(line[30:38].strip()),
                        'y': float(line[38:46].strip()),
                        'z': float(line[46:54].strip()),
                        'occupancy': float(line[54:60].strip() or '1.00'),
                        'temp_factor': float(line[60:66].strip() or '0.00'),
                        'element': line[76:78].strip(),
                        'charge': line[78:80].strip()
                    }
                    atoms.append(atom)
                    serial_counter += 1
                except ValueError:
                    continue  # Skip lines with parsing issues
    return atoms

def write_pdb_atoms(atoms, output_file_path):
    # --- Minimal change: renumber atoms sequentially ---
    for i, atom in enumerate(atoms, start=1):
        atom['serial'] = i
    # ------------------------------------------------------
    with open(output_file_path, 'w') as f:
        for atom in atoms:
            f.write(format_pdb_line(atom))
        f.write("END\n")

def format_pdb_line(atom):
    line = "{:<6}{:>5} {:<4}{:1}{:<3} {:1}{:>4}{:<1}   ".format(
        atom['record_name'], atom['serial'],
        atom['atom_name'].ljust(4), atom['alt_loc'],
        atom['res_name'], atom['chain_id'],
        atom['res_seq'], atom['i_code']
    )
    line += "{:>8.3f}{:>8.3f}{:>8.3f}".format(
        atom['x'], atom['y'], atom['z']
    )
    line += "{:>6.2f}{:>6.2f}          ".format(
        atom['occupancy'], atom['temp_factor']
    )
    line += "{:<2}{:<2}\n".format(atom['element'], atom['charge'])
    return line

def setup_logging(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

if __name__ == "__main__":
    main()

