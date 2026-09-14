# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 3
"""
Usage (from terminal):
    python Exercise_3.py structure.pdb
    python Exercise_3.py structure.pdb --cutoff 3.2
"""

import argparse
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

polar_atoms = ('O', 'N', 'S')

def get_label(residue):
    seg = residue.get_parent().id
    num = residue.id[1]
    ins = residue.id[2].strip()
    return f"{residue.get_resname()} {seg}{num}{ins}"


def main():
    argp = argparse.ArgumentParser(
        prog='Exercise_3',
        description='List possible hydrogen bonds between polar atoms'
    )

    argp.add_argument(
        '--cutoff',
        dest='cutoff',
        type=float,
        default=3.5,
        help='Distance criterium (dist < 3.5 Å)'
    )

    argp.add_argument(
        'input_pdb',
        help='PDB file'
    )

    opts = argp.parse_args()

    struct_name = os.path.splitext(os.path.basename(opts.input_pdb))[0]

    reader = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = reader.get_structure(struct_name, opts.input_pdb)

    polar_list = [a for a in structure.get_atoms() if a.element in polar_atoms]

    searcher = NeighborSearch(polar_list)

    bonds = []
    for a1, a2 in searcher.search_all(opts.cutoff):
        r1 = a1.get_parent()
        r2 = a2.get_parent()
        if r1 is r2:
            continue
        d = a1 - a2
        bonds.append((r1, a1, r2, a2, d))

    bonds.sort(key=lambda h: (h[0].get_parent().id, h[0].id[1], h[2].get_parent().id, h[2].id[1]))

    print(f"Possible hydrogen bonds (polar atom-atom distance < {opts.cutoff} A)")
    print("-" * 70)
    for r1, a1, r2, a2, d in bonds:
        print(f"{get_label(r1):>12}.{a1.get_name():<4} -- "
              f"{get_label(r2):<12}.{a2.get_name():<4}  {d:6.2f} A")

    print(f"\nTotal possible hydrogen bonds found: {len(bonds)}")


if __name__ == '__main__':
    main()