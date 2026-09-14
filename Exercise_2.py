# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 2
"""
Usage (from terminal):
    python Exercise_2.py structure.pdb 35
    python Exercise_2.py --chain A structure.pdb 35
"""

import argparse
import os
import sys

from Bio.PDB.PDBParser import PDBParser


def main():
    argp = argparse.ArgumentParser(
        prog='Exercise_2',
        description='List all atoms (name and coordinates) for a given residue number')

    argp.add_argument(
        '--chain',
        dest='seg_id',
        default=None,
        help='Chain id (optional). If not given, all chains are searched')

    argp.add_argument(
        'input_pdb',
        help='Input PDB file')

    argp.add_argument(
        'res_num',
        type=int,
        help='Residue number')

    opts = argp.parse_args()

    struct_name = os.path.splitext(os.path.basename(opts.input_pdb))[0]

    reader = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = reader.get_structure(struct_name, opts.input_pdb)

    model0 = structure[0]

    hits = []
    for seg in model0:
        if opts.seg_id is not None and seg.id != opts.seg_id:
            continue
        for res in seg:
            if res.id[1] == opts.res_num:
                hits.append(res)

    if not hits:
        sys.exit(f"No residue with number {opts.res_num} found"
                  f"{' in chain ' + opts.seg_id if opts.seg_id else ''}.")

    for res in hits:
        seg_id = res.get_parent().id
        print(f"\nResidue: {res.get_resname()} {seg_id}{res.id[1]}")
        print("-" * 50)
        print(f"{'Atom':<6}{'X':>10}{'Y':>10}{'Z':>10}")
        atom_list = sorted(res.get_atoms(), key=lambda a: a.get_serial_number())
        for a in atom_list:
            x, y, z = a.get_coord()
            print(f"{a.get_name():<6}{x:>10.3f}{y:>10.3f}{z:>10.3f}")


if __name__ == '__main__':
    main()