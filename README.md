# Restore Protein Structure Pipeline

## Overview

This repository contains a Python-based pipeline for restoring and refining protein structures in **PDB format**. The pipeline follows a three-step process:

1. **Add missing heavy atoms** using [PDBFixer](https://github.com/openmm/pdbfixer).
2. **Optimize side-chain conformations** using [SCWRL4](http://dunbrack.fccc.edu/scwrl4/).
3. **Merge atomic coordinates** from the original, fixed, and optimized structures while preserving **beta-factor information**.

This tool is designed for **structural bioinformatics** applications, enabling researchers to reconstruct high-quality protein models from incomplete structures.

---

## Features

- **Automated heavy atom reconstruction** for incomplete PDB files.
- **Parallel processing** for high throughput.
- **Side-chain optimization** with SCWRL4.
- **Preservation of atomic beta-factors** during reconstruction.
- **Detailed logging and statistical reports** for quality assessment.

---

## Dependencies

### Required Python Packages

- **Python 3.x**
- **PDBFixer (v1.10.0)**
- **OpenMM (v8.2.0)**
- **NumPy (v2.2.0rc1)**
- **tqdm (v4.67.1)**

### External Software

- **SCWRL4** – A required executable for side-chain optimization.
  - Ensure `Scwrl4` is installed and available in your system's `PATH`.
  - Download SCWRL4 from [here](http://dunbrack.fccc.edu/scwrl4/).

---

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/restore-protein-structure.git
cd restore-protein-structure
```

### 2. Create a Virtual Enviroment (Optional but Recommended)

```bash
python3 -m venv venv
source venv/bin/activate  # On macOS/Linux
venv\Scripts\activate     # On Windows
```

### 3. Install Python Dependencies
```bash
pip install -r requirements.txt
```
or install dependencies manually:

```bash
  pip install pdbfixer openmm numpy tqdm
```

### 4. Verify SCWRL4 Installation

```bash
Scwrl4 -h
```

If the command is not found, add the SCWRL4 directory to your PATH.

## Usage

### Running the Pipeline

```bash
./restore_prot_struct.py --input /path/to/input_pdbs --output /path/to/output_directory
```

### Command-Line Arguments

| Argument | Description |
|------------|--------------------------------------------------|
| `--input` | Path to the directory containing input PDB files | 
| `--output` | Path to the directory for output PDB files |

The script processes all `.pdb` files in the input directory.

## Output Files

For each input PDB file (e.g., `example.pdb`), the pipeline generates:

| Output File | Description |
|----------------------|-------------|
| `example_A.pdb` | PDB after adding missing heavy atoms with PDBFixer. |
| `example_B.pdb` | PDB after side-chain optimization with SCWRL4. |
| `example_C.pdb` | Final merged PDB combining the original, fixed, and optimized structures. |

### Additional Files:
- `fix_pdbs.log` – Logs the processing details.
- `stats.dat` – Statistical summary of missing atoms per file.

## Logging and Statistics

- Logging:
Processing details and errors are logged in `fix_pdbs.log` within the output directory.

- Statistics:
The script computes:
- Total files processed
- Successful/failed processing counts
- Average missing atoms per file
- Standard deviation of missing atoms

These statistics are printed to the console and saved to `stats.dat`.

## Example Workflow

### Step 1: Prepare Input Files
Place all `.pdb` files in a designated input directory:

```
/home/user/input_pdbs/
├── protein1.pdb
├── protein2.pdb
├── protein3.pdb
```

### Step 2: Run the Script

```bash
./restore_prot_struct.py --input /home/user/input_pdbs --output /home/user/output_pdbs
```

### Step 3: Check Output Files

```
/home/user/output_pdbs/
├── protein1_A.pdb
├── protein1_B.pdb
├── protein1_C.pdb
├── fix_pdbs.log
├── stats.dat
```

## Troubleshooting

### 1. SCWRL4 Not Found

**Error:**
```
Error: SCWRL4 executable 'Scwrl4' not found in PATH.
```

**Solution:**
- Verify SCWRL4 is installed:

```bash
Scwrl4 -h
```

- If not found, add its location to `PATH`:

```bash
export PATH="/path/to/Scwrl4:$PATH"
```

### 2. Missing Dependencies

**Error:**
```
ModuleNotFoundError: No module named 'pdbfixer'
```

**Solution:**
Ensure dependencies are installed:

```bash
pip install pdbfixer openmm numpy tqdm
```

### 3. No Output Files Generated

- Verify `.pdb` files exist in the input directory.
- Check `fix_pdbs.log` for detailed error messages.

## License

This project is distributed under the MIT License.


## Acknowledgments

PDBFixer
- PDBFixer(https://github.com/openmm/pdbfixer) – Used for fixing PDB files.
  
OpenMM
- OpenMM(https://github.com/openmm/openmm) – Molecular simulation toolkit.
  
SCWRL4
- SCWRL4(http://dunbrack.fccc.edu/scwrl4/) – Side-chain optimization tool.
