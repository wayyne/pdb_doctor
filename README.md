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

