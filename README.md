# Restore Protein Structure Pipeline

This repository contains a Python pipeline designed for restoring and refining protein structures in PDB format. The pipeline performs the following operations on input PDB files:

**Addition of missing heavy atoms** using 
P
D
B
F
i
x
e
r
PDBFixer(https://github.com/openmm/pdbfixer)
**Side-chain optimization** with 
S
C
W
R
L
4
SCWRL4(http://dunbrack.fccc.edu/scwrl4/)
**Merging of coordinates** from the original, fixed, and optimized structures while transferring beta-factor information
The code is suitable for academic research purposes where accurate protein structure restoration is required.

---

## Overview

The main script (`restore_prot_struct.py`) processes all PDB files found in a given input directory. For each file, it:

**Detects and adds missing heavy atoms** with PDBFixer.
**Optimizes side-chain conformations** by invoking the SCWRL4 executable.
**Merges coordinates** to produce a final output file that preserves beta-factor information from the original structure.
The pipeline executes these steps in parallel for improved performance and provides detailed logging and statistical summaries upon completion.

---

## Dependencies

The pipeline depends on the following packages and external tools:

**Python 3.x**
**PDBFixer (v1.10.0)**
**OpenMM (v8.2.0)**
**NumPy (v2.2.0rc1)**
**tqdm (v4.67.1)**
**SCWRL4** - This executable must be available in your system's PATH. (Installation instructions for SCWRL4 can be found on its 
o
f
f
i
c
i
a
l
w
e
b
s
i
t
e
officialwebsite(http://dunbrack.fccc.edu/scwrl4/).)
Additional standard Python libraries used include `os`, `sys`, `argparse`, `subprocess`, `collections`, `shutil`, `concurrent.futures`, `logging`, and `statistics`.

---

## Installation

**Clone the repository:**
``` git clone https://github.com/yourusername/restore-protein-structure.git cd restore-protein-structure ```
**Set up a virtual environment (optional but recommended):**
``` python3 -m venv venv source venv/bin/activate ```
**Install the required Python packages:**
``` pip install -r requirements.txt ```

*Note:* If a `requirements.txt` file is not provided, install the dependencies manually:

``` pip install pdbfixer openmm numpy tqdm ```
**Ensure SCWRL4 is installed and available in your PATH:**
Verify by running:

``` Scwrl4 -h ```
---

## Usage

Run the pipeline by executing the main script with the required input and output directory arguments. For example:

``` ./restore_prot_struct.py --input /path/to/input_pdbs --output /path/to/output_directory ```

**Command-Line Arguments**

`--input`: Path to the directory containing input PDB files.
`--output`: Path to the directory where the processed PDB files and logs will be saved.
The script processes all files with the `.pdb` extension found in the input directory.

---

## Output Files

For each input PDB file (e.g., `example.pdb`), the pipeline produces three output files in the specified output directory:

**example_A.pdb**: Result after adding missing heavy atoms with PDBFixer.
**example_B.pdb**: Result after side-chain optimization via SCWRL4.
**example_C.pdb**: Final merged PDB combining coordinates and beta-factor information from the original, fixed, and optimized structures.
Additionally, a log file (`fix_pdbs.log`) and a statistics file (`stats.dat`) are created in the output directory summarizing the processing results.

---

## Logging and Statistics

**Logging:**
The pipeline logs detailed processing information and error messages to `fix_pdbs.log` in the output directory.
**Statistics:**
After processing, the script computes and logs:
The total number of files processed.
Counts of successful and failed processing attempts.
The average and standard deviation of missing heavy atoms per file.
These statistics are printed to the console and saved to `stats.dat`.

---

## Contributing

Contributions are welcome. If you have suggestions or improvements, please submit a pull request or open an issue on the GitHub repository.

---

## License

This project is distributed under the 
M
I
T
L
i
c
e
n
s
e
MITLicense(LICENSE).

---

## Acknowledgements

**PDBFixer** and **OpenMM** for providing essential tools for protein structure processing.
**SCWRL4** for its robust side-chain prediction capabilities.
The developers and maintainers of the underlying libraries and tools.
---

## Contact

For further questions or comments regarding this pipeline, please contact
Y
o
u
r
N
a
m
e
YourName at 
y
o
u
r
.
e
m
a
i
l
@
e
x
a
m
p
l
e
.
c
o
m
your.email@example.com.
