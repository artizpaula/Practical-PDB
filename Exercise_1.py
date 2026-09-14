# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 1
"""
Usage (from terminal):
    python Exercise_1.py structure.pdb 5.0
    python Exercise_1.py --dist 5.0 structure.pdb
"""

import argparse
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch


def get_label(residue):
    seg = residue.get_parent().id
    num = residue.id[1]
    ins = residue.id[2].strip()
    return f"{residue.get_resname()} {seg}{num}{ins}"


def main():
    argp = argparse.ArgumentParser(
        prog='Exercise_1',
        description='List pairs of residues whose CA atoms are closer than a given distance')

    argp.add_argument(
        '--dist',
        dest='cutoff',
        type=float,
        default=5.0,
        help='Distance threshold in Angstroms (default: 5.0)')

    argp.add_argument(
        'input_pdb',
        help='Input PDB file')

    opts = argp.parse_args()

    struct_name = os.path.splitext(os.path.basename(opts.input_pdb))[0]

    reader = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = reader.get_structure(struct_name, opts.input_pdb)

    ca_list = [a for a in structure.get_atoms() if a.id == 'CA']

    searcher = NeighborSearch(ca_list)

    found = []
    for a1, a2 in searcher.search_all(opts.cutoff):
        r1 = a1.get_parent()
        r2 = a2.get_parent()
        if r1 is r2:
            continue
        d = a1 - a2
        found.append((r1, r2, d))

    found.sort(key=lambda p: (p[0].get_parent().id, p[0].id[1], p[1].get_parent().id, p[1].id[1]))

    print(f"Pairs of residues with CA-CA distance < {opts.cutoff} A")
    print("-" * 60)
    for r1, r2, d in found:
        print(f"{get_label(r1):>12} -- {get_label(r2):<12}  {d:6.2f} A")

    print(f"\nTotal pairs found: {len(found)}")


if __name__ == '__main__':
    main()