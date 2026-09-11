# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python
#
""" Exercise 2
Generate a list of all atoms for a given residue number.

Parameters: PDB file name, Residue number (Including Chain if applicable)

Usage (from terminal):
    python Exercise_2.py structure.pdb 35
    python Exercise_2.py --chain A structure.pdb 35
"""

import argparse
import os
import sys

from Bio.PDB.PDBParser import PDBParser


def main():
    parser = argparse.ArgumentParser(
        prog='Exercise_2',
        description='List all atoms (name and coordinates) for a given residue number'
    )

    parser.add_argument(
        '--chain',
        dest='chain',
        default=None,
        help='Chain id (optional). If not given, all chains are searched'
    )

    parser.add_argument(
        'pdb_file',
        help='Input PDB file'
    )

    parser.add_argument(
        'resnum',
        type=int,
        help='Residue number'
    )

    args = parser.parse_args()

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]

    pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)
    st = pdb_parser.get_structure(pdb_id, args.pdb_file)

    model = st[0]

    matched_residues = []
    for chain in model:
        if args.chain is not None and chain.id != args.chain:
            continue
        for res in chain:
            if res.id[1] == args.resnum:
                matched_residues.append(res)

    if not matched_residues:
        sys.exit(f"No residue with number {args.resnum} found"
                  f"{' in chain ' + args.chain if args.chain else ''}.")

    for res in matched_residues:
        chain_id = res.get_parent().id
        print(f"\nResidue: {res.get_resname()} {chain_id}{res.id[1]}")
        print("-" * 50)
        print(f"{'Atom':<6}{'X':>10}{'Y':>10}{'Z':>10}")
        # Sort atoms by serial number for a consistent, easy-to-read order
        atoms = sorted(res.get_atoms(), key=lambda at: at.get_serial_number())
        for at in atoms:
            x, y, z = at.get_coord()
            print(f"{at.get_name():<6}{x:>10.3f}{y:>10.3f}{z:>10.3f}")


if __name__ == '__main__':
    main()