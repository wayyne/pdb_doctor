# PDB Doctor

## Overview

**PDB Doctor** is a powerful, **Unix-driven** pipeline designed for **protein structure refinement** and **completion**. It automates the extraction, folding, and restoration of protein structures from **metadata-scrubbed PDB files**, utilizing **PDBFixer, SCWRL4, and ESM3** for **atom reconstruction, backbone restoration, and sequence-based folding**.

This tool is particularly useful for:
- **Extracting accurate amino acid sequences** from incomplete PDB files.
- **Folding sequences using ESM3** and aligning them to experimental data.
- **Restoring missing residues** while ensuring atomic correctness.
- **Completing PDB structures** with side-chain optimization via SCWRL4.
- **Sanity-checking amide bonds** and atom consistency.

## Features
- **Automated Extraction:** Retrieves protein sequences from metadata-scrubbed PDBs.
- **Folding & Alignment:** Uses **ESM3** to predict missing regions and align them globally.
- **Local Stitching:** Inserts missing residues with **backbone preservation**.
- **Heavy Atom Completion:** Fixes missing **heavy atoms** using **PDBFixer**.
- **Side-Chain Optimization:** Restores **side-chain conformations** using **SCWRL4**.
- **Fully Scripted Workflow:** A single **Bash script (`triage.sh`)** automates the pipeline.

---

## Installation

PDB Doctor requires a **Python 3.10+ environment**. You can set up dependencies using **either Conda or Python venv** via the provided script.

### 1. Clone the Repository

```bash
git clone https://github.com/wayyne/pdb_doctor.git
cd pdb_doctor
```

### 2. Install Dependencies

Use the **setup script** to configure either **Conda** or **venv**:

#### **Option A: Conda Setup**
```bash
cd env
bash setup_env.sh conda
```

#### **Option B: Python venv Setup**
```bash
cd env
bash setup_env.sh venv
```

This will:
- Create a Python environment (`pdb_doctor`).
- Install dependencies from `pdb_doctor_env.txt`.
- Set up **PDBFixer** (from OpenMM) and **SCWRL4** (if available).

### 3. Ensure SCWRL4 is Installed
SCWRL4 is required for **side-chain optimization**. Make sure it is in your `PATH`:
```bash
which Scwrl4
```
If not installed, download it from:
[http://dunbrack.fccc.edu/scwrl4/](http://dunbrack.fccc.edu/lab/scwrl/)

---

## Usage

### **1. Run the Full Workflow (assuming python venv)**
The **`triage.sh`** script automates the complete process:

```bash
source <pdb_doctor_root>/env/pdb_doctor/bin/activate
bash triage.sh
```

### **2. Workflow Breakdown**
Below is a breakdown of how the pipeline functions.

#### **Step 1: Extract Sequences**
Extracts **full-length amino acid sequences** from metadata-scrubbed PDBs.
```bash
python extract_sequence.py <input_directory>
```
**Outputs:**
- `complete.tsv` – Fully resolved sequences.
- `incomplete.tsv` – Sequences with **missing residues**.

---

#### **Step 2: Fold and Align Structures**
Folds missing segments using **ESM3**, aligns the predicted structure, and **stitches missing regions** into the experimental model.
```bash
python fold_and_transfer.py <partial_pdb> "<sequence>" <output_pdb>
```
**Outputs:**
- `fitted_folded.pdb` – Folded model, globally aligned.
- `output_pdb` – Merged **chimeric structure** with restored residues.

---

#### **Step 3: Complete Missing Heavy Atoms**
Uses **PDBFixer** to **restore heavy atoms** and **SCWRL4** to refine side-chain placements.
```bash
python complete_pdb.py --mode partial --input input --output output
```
**Modes:**
- `--mode partial` → Restores missing **heavy atoms**.

---

## Example Workflow

Given a directory containing **scrubbed PDB files**, PDB Doctor can:
1. Extract sequences.
2. Fold missing regions using **ESM3**.
3. Align folded models with **Kabsch-based superposition**.
4. Locally **stitch missing segments**.
5. Restore **side-chains** using **SCWRL4**.

Example Run:
```bash
bash triage.sh
```

After execution:
- `output/` will contain **corrected PDB structures**.
- `extraction.log` will provide **sequence recovery details**.

---

## Dependencies
- **Python 3.10+**
- **NumPy**
- **Torch** (for ESM3 folding)
- **pdbfixer** (for atom reconstruction)
- **SCWRL4** (for side-chain refinement)
- **Requests** (for fetching PDBs)

---

## Citation
If you use PDB Doctor in your research, please cite:
- **SCWRL4:** Krivov et al., *Proteins: Structure, Function, and Bioinformatics (2009).*
- **PDBFixer:** OpenMM PDBFixer repository.

---

## Author
Guy **Wayyne** Dayhoff II

---

## License
MIT License.  
Feel free to contribute or report issues!

