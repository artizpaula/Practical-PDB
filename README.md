# Practical PDB

Practical 1 exercises for **Introduction to protein structure manipulation in Python**, using [Biopython](https://biopython.org/)'s `Bio.PDB` module to parse, inspect and analyze protein structures in PDB format.

**Authors:** Alicia Mañas, Lídia Sanchez, Paula Artiz

## Repository contents

| File | Description |
|---|---|
| `Exercise_1.py` | Lists pairs of residues whose CA atoms are closer than a given distance |
| `Exercise_2.py` | Lists all atoms (name and coordinates) for a given residue number |
| `Exercise_3.py` | Lists possible hydrogen bonds between polar atoms (O, N, S) |
| `Exercise_4.py` | Lists all CA atoms of a given residue type (1- or 3-letter code) |
| `Exercise_5.py` | Lists possible peptide bonds (backbone C–N contacts) |
| `Exercise_6.py` | Lists possible disulphide bonds (Cys SG–SG contacts) |
| `Exercise_7.py` | Prints distances between all atom pairs of two given residues |
| `Exercise_8.py` | Optional: download a structure by PDB id instead of using a local file |
| `Exercise_*.txt` | Example command-line inputs and outputs for each exercise |
| `1UBQ.pdb`, `4HHB.pdb` | Example structures used to test the scripts |

## Usage

Each script uses `argparse`; run with `-h` to see all options. General pattern:

```bash
python Exercise_1.py <input_pdb> <distance>
python Exercise_2.py [--chain CHAIN] <input_pdb> <res_num>
python Exercise_3.py [--cutoff CUTOFF] <input_pdb>
python Exercise_4.py <input_pdb> <res_type>
python Exercise_5.py [--cutoff CUTOFF] <input_pdb>
python Exercise_6.py [--cutoff CUTOFF] <input_pdb>
```

### Examples

```bash
# Residue pairs with CA atoms closer than 5 A
python Exercise_1.py 1UBQ.pdb 5

# All atoms of residue 40
python Exercise_2.py 4HHB.pdb 40

# Possible hydrogen bonds (default cutoff 3.5 A)
python Exercise_3.py 1UBQ.pdb

# All ARG residues (CA atoms)
python Exercise_4.py 4HHB.pdb ARG
```
