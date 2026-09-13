#!/usr/bin/env python

""" Exercise 5
Ordinary peptide bonds are made between atom C of one residue and atom N
of the following. Usual distance should be below 2 A.
Find pairs of C-N atoms from different residues that are closer than
the cut-off distance (default 2.0 A).

Parameters: PDB file name. Optional: cut-off distance (defaults to 2.0)

Usage (from terminal):
    python Exercise_5.py structure.pdb
    python Exercise_5.py structure.pdb --cutoff 2.2
"""

import argparse
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch


def residue_id(res):
    chain_id = res.get_parent().id
    resnum = res.id[1]
    icode = res.id[2].strip()
    return f"{res.get_resname()} {chain_id}{resnum}{icode}"


def main():
    parser = argparse.ArgumentParser(
        prog='Exercise_5',
        description='List possible peptide bonds (C-N contacts) between residues'
    )

    parser.add_argument(
        '--cutoff',
        dest='cutoff',
        type=float,
        default=2.0,
        help='Distance criterium (dist < 2.0 A)'
    )

    parser.add_argument(
        'pdb_file',
        help='PDB file'  # 1UBQ.pdb or 4HHB.pdb
    )

    args = parser.parse_args()

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]

    pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)
    st = pdb_parser.get_structure(pdb_id, args.pdb_file)

    # Select only backbone C and N atoms (identified by atom name, not element)
    cn_atoms = [at for at in st.get_atoms() if at.get_name() in ('C', 'N')]

    nbsearch = NeighborSearch(cn_atoms)

    pbonds = []
    for at1, at2 in nbsearch.search_all(args.cutoff):
        # Only interested in C...N pairs (one of each), not C-C or N-N
        names = {at1.get_name(), at2.get_name()}
        if names != {'C', 'N'}:
            continue

        # Make sure at1 is always the C and at2 is always the N, for readability
        if at1.get_name() != 'C':
            at1, at2 = at2, at1

        res1 = at1.get_parent()
        res2 = at2.get_parent()

        # A peptide bond only makes sense between atoms of different residues
        if res1 is res2:
            continue

        dist = at1 - at2
        pbonds.append((res1, at1, res2, at2, dist))

    # Sort by residue number / chain of the first atom involved
    pbonds.sort(key=lambda h: (h[0].get_parent().id, h[0].id[1], h[2].get_parent().id, h[2].id[1]))

    print(f"Possible peptide bonds (C-N distance < {args.cutoff} A)")
    print("-" * 70)
    for res1, at1, res2, at2, dist in pbonds:
        print(f"{residue_id(res1):>12}.{at1.get_name():<4} -- "
              f"{residue_id(res2):<12}.{at2.get_name():<4}  {dist:6.2f} A")

    print(f"\nTotal possible peptide bonds found: {len(pbonds)}")


if __name__ == '__main__':
    main()