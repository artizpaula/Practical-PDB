# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 4
"""
Usage (from terminal):
    python Exercise_4.py structure.pdb ARG
    python Exercise_4.py structure.pdb R
"""

import argparse
import os
import sys

from Bio.PDB.PDBParser import PDBParser

AA3_TO_AA1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}
AA1_TO_AA3 = {one: three for three, one in AA3_TO_AA1.items()}


def to_three_letter(code):
    code = code.strip().upper()
    if len(code) == 1:
        if code not in AA1_TO_AA3:
            sys.exit(f"Unknown one-letter residue code: {code}")
        return AA1_TO_AA3[code]
    elif len(code) == 3:
        if code not in AA3_TO_AA1:
            sys.exit(f"Unknown three-letter residue code: {code}")
        return code
    else:
        sys.exit(f"Residue type must be a 1- or 3-letter code, got: {code}")


def main():
    argp = argparse.ArgumentParser(
        prog='Exercise_4',
        description='List all CA atoms (with coordinates) of a given residue type')

    argp.add_argument(
        'input_pdb',
        help='Input PDB file')

    argp.add_argument(
        'res_type',
        help='Residue type, either one-letter (e.g. R) or three-letter (e.g. ARG) code')

    opts = argp.parse_args()

    target_resname = to_three_letter(opts.res_type)

    struct_name = os.path.splitext(os.path.basename(opts.input_pdb))[0]

    reader = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = reader.get_structure(struct_name, opts.input_pdb)

    matches = []
    for a in structure.get_atoms():
        if a.id == 'CA' and a.get_parent().get_resname() == target_resname:
            matches.append(a)

    matches.sort(key=lambda a: (a.get_parent().get_parent().id, a.get_parent().id[1]))

    print(f"CA atoms of residue type {target_resname} ({AA3_TO_AA1[target_resname]})")
    print("-" * 60)
    print(f"{'Chain':<6}{'ResNum':<8}{'X':>10}{'Y':>10}{'Z':>10}")
    for a in matches:
        res = a.get_parent()
        seg_id = res.get_parent().id
        x, y, z = a.get_coord()
        print(f"{seg_id:<6}{res.id[1]:<8}{x:>10.3f}{y:>10.3f}{z:>10.3f}")

    print(f"\nTotal {target_resname} residues found: {len(matches)}")


if __name__ == '__main__':
    main()