# 🩺 PDB Doctor

## Overview

**PDB Doctor** is a powerful, **Unix-driven** pipeline designed for **protein structure refinement** and **completion**.  
It automates the **extraction**, **folding**, **transplantation**, and **restoration** of protein structures from **metadata-scrubbed PDB files**, utilizing:

- **PDBFixer** for atom reconstruction  
- **SCWRL4** for side-chain refinement  
- **ESM3** for sequence-based folding and structure completion. (ESM Forge token required: place in pdb_doctor/env/forge_token.txt)

This tool is particularly useful for:
- Extracting accurate amino acid sequences from incomplete PDB files.  
- Folding sequences using ESM3 and aligning them to experimental data.  
- Restoring missing residues while ensuring atomic correctness.  
- Completing PDB structures with side-chain optimization via SCWRL4.  
- Sanity-checking amide bonds and atom consistency.

---

## 🧰 Core Tools

| Command | Role | Description | When to use |
|---------|------|-------------|-------------|
| \`pdbdr-triage\` | **Triage Nurse** | Batch intake of PDBs. Classifies complete vs incomplete, routes incomplete through folding or transplant, runs completion, and contact labeling. | Full directory processing |
| \`pdbdr-transplant\` | **Transplant Surgeon** | Fills missing residues using a donor PDB (via PDB ID) and grafts into the partial structure. | If a suitable donor exists |
| \`pdbdr-bioprint\` | **Biofabrication** | Folds a structure de novo from a sequence and grafts missing regions into the partial structure. | When no donor PDB exists |
| \`pdbdr-complete-pdb\` | **Finishing Step** | Restores missing heavy atoms and refines side chains using PDBFixer and SCWRL4. | Always after grafting |
| \`pdbdr-contact-trace\` | **Contact Tracer** | Reattaches ligands (HETATMs) and labels contacting residues. | Post-completion |
| \`pdbdr-extract-seq\` | **Sequence Extraction** | Extracts sequences from metadata-scrubbed PDB files. | Pre-triage or standalone |

---

## ✨ Key Features

- Automated sequence extraction  
- De novo or donor-based structure completion  
- Global alignment (Kabsch) and local stitching  
- Heavy atom reconstruction via PDBFixer  
- **SCWRL4-based side-chain refinement (required)**  
- Ligand contact labeling and annotation  
- Fully scripted, cluster-friendly workflow

---

## 📦 Installation

PDB Doctor requires **Python 3.10+**.

### 1. Clone the Repository
```
git clone https://github.com/wayyne/pdb_doctor.git
cd pdb_doctor
```

---

### 2. Choose an Installation Route

PDB Doctor can be set up in **two main ways**, depending on how you want to use it:

#### **Option A: Recommended — Dedicated venv Environment**

```
cd env
bash setup_env.sh venv
```

✅ Recommended because:
- Creates a clean, isolated environment
- Avoids Conda/Torch conflicts
- All \`pdbdr-*\` command-line tools are automatically added to your PATH
- You can run tools directly as binaries (e.g. \`pdbdr-triage\`)

⚠️ **Important:**  
If you use Conda for other work, make sure it’s fully deactivated before creating or activating a Python venv:
```
conda deactivate
```

---

#### **Option B: Conda Environment (Not Recommended)**

```
cd env
bash setup_env.sh conda
```

⚠️ Known issues:  
- Conda can sometimes conflict with ESM3 or Torch builds.  
- You may need to resolve dependencies manually.

---

#### **Option C: Install into an Existing Environment or System-Wide**

```
bash src/install_pdbdr.sh
```

This:
- Installs PDB Doctor as a **Python package** in the currently active environment  
- Or system-wide if no virtual environment is active  
- **Does not** create a new environment

✅ This allows you to:
- \`import pdbdr\` from within your own Python projects or notebooks  
- Use \`pdbdr\` functions programmatically as a library  
- Integrate the tools into shared lab or cluster environments

📝 If you want the **CLI binaries** (e.g., \`pdbdr-transplant\`, \`pdbdr-triage\`) to work from anywhere, make sure the Python environment’s \`bin\` directory is on your \`PATH\`.

---

### 3. **Install SCWRL4 (Required)**

SCWRL4 is required for side-chain optimization and must be available in your PATH:
```
which Scwrl4
```

If not installed, download it from:
[http://dunbrack.fccc.edu/scwrl4/](http://dunbrack.fccc.edu/lab/scwrl/)

Place the binary somewhere in your PATH (e.g., \`~/bin\` or \`/usr/local/bin\`).

---

## 🧪 Usage

### 1. Run the Full Workflow
```
source env/pdb_doctor/bin/activate
pdbdr-triage
```

This command:
- Extracts sequences
- Detects incomplete structures
- Runs transplant or bioprint
- Completes missing atoms
- Labels ligand contacts
- Organizes triaged outputs

---

### 2. Individual Steps

#### Extract Sequences
```
pdbdr-extract-seq .
```

#### Transplant Missing Segments
```
pdbdr-transplant partial.pdb 1ABC output.pdb --chain A
```

#### Bioprint from Sequence
```
pdbdr-bioprint partial.pdb "MSTK..." output.pdb
```

#### Complete Missing Heavy Atoms and Side Chains
```
pdbdr-complete-pdb --mode partial --input input_dir --output output_dir
```

#### Label Ligand Contacts
```
pdbdr-contact-trace input_dir output_dir
```

---

## 🧭 Typical Workflow

1. **Triage** — classify PDBs as complete or incomplete  
2. **Transplant / Bioprint** — restore missing residues  
3. **Complete** — rebuild heavy atoms and refine side chains (**requires SCWRL4**)  
4. **Contact Trace** — label ligand contacts  
5. **Export** — output fully completed chimeric structures

---

## 🧬 Dependencies

- Python 3.10+  
- NumPy  
- Torch (for ESM3 folding)  
- PDBFixer (atom reconstruction)  
- **SCWRL4 (required)** — for side-chain refinement  
- Requests (fetching PDBs)

---

## 📜 Citation

If you use **PDB Doctor** in your research, please cite:

- Krivov et al., *Proteins: Structure, Function, and Bioinformatics (2009)* — SCWRL4  
- OpenMM PDBFixer repository  
- ESM3 model (Meta AI)

---

## 👨‍⚕️ Author

Guy **Wayyne** Dayhoff II

---

## 🪪 License

MIT License — Contributions welcome!

