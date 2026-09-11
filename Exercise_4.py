# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python
#
""" Exercise 4
Generate a list of all CA atoms of given residue type with coordinates.

Parameters: PDB file name, residue type.
Optional: accept residue codes in one- or three-letter formats automatically

Usage (from terminal):
    python Exercise_4.py structure.pdb ARG
    python Exercise_4.py structure.pdb R
"""

import argparse
import os
import sys

from Bio.PDB.PDBParser import PDBParser

# Standard 20 amino acids: three-letter <-> one-letter code
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}
ONE_TO_THREE = {one: three for three, one in THREE_TO_ONE.items()}


def normalize_resname(code):
    """Accept either a one-letter or a three-letter residue code and
    return the three-letter code used in PDB files."""
    code = code.strip().upper()
    if len(code) == 1:
        if code not in ONE_TO_THREE:
            sys.exit(f"Unknown one-letter residue code: {code}")
        return ONE_TO_THREE[code]
    elif len(code) == 3:
        if code not in THREE_TO_ONE:
            sys.exit(f"Unknown three-letter residue code: {code}")
        return code
    else:
        sys.exit(f"Residue type must be a 1- or 3-letter code, got: {code}")


def main():
    parser = argparse.ArgumentParser(
        prog='Exercise_4',
        description='List all CA atoms (with coordinates) of a given residue type'
    )

    parser.add_argument(
        'pdb_file',
        help='Input PDB file'
    )

    parser.add_argument(
        'restype',
        help='Residue type, either one-letter (e.g. R) or three-letter (e.g. ARG) code'
    )

    args = parser.parse_args()

    resname = normalize_resname(args.restype)

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]

    pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)
    st = pdb_parser.get_structure(pdb_id, args.pdb_file)

    selected = []
    for at in st.get_atoms():
        if at.id == 'CA' and at.get_parent().get_resname() == resname:
            selected.append(at)

    # Sort by chain, then residue number
    selected.sort(key=lambda at: (at.get_parent().get_parent().id, at.get_parent().id[1]))

    print(f"CA atoms of residue type {resname} ({THREE_TO_ONE[resname]})")
    print("-" * 60)
    print(f"{'Chain':<6}{'ResNum':<8}{'X':>10}{'Y':>10}{'Z':>10}")
    for at in selected:
        res = at.get_parent()
        chain_id = res.get_parent().id
        x, y, z = at.get_coord()
        print(f"{chain_id:<6}{res.id[1]:<8}{x:>10.3f}{y:>10.3f}{z:>10.3f}")

    print(f"\nTotal {resname} residues found: {len(selected)}")


if __name__ == '__main__':
    main()