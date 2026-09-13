#!/usr/bin/env python

""" Exercise 6
Disulphide bonds are formed between S atoms (SG) of Cys residues when they
are at the appropriate distance (around 1.9 A). A wider cut-off is allowed
by default to account for structural variability.

Parameters: PDB file name. Optional: cut-off distance (defaults to 2.5)

Usage (from terminal):
    python Exercise_6.py structure.pdb
    python Exercise_6.py structure.pdb --cutoff 2.2
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
        prog='Exercise_6',
        description='List possible disulphide bonds (Cys SG-SG contacts)'
    )

    parser.add_argument(
        '--cutoff',
        dest='cutoff',
        type=float,
        default=1.9,
        help='Distance criterium (dist < 1.9 A)'
    )

    parser.add_argument(
        'pdb_file',
        help='PDB file'  # 1UBQ.pdb or 4HHB.pdb
    )

    args = parser.parse_args()

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]

    pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)
    st = pdb_parser.get_structure(pdb_id, args.pdb_file)

    # Select only the side-chain sulfur (SG) atoms of Cys residues
    sg_atoms = [
        at for at in st.get_atoms()
        if at.get_parent().get_resname() == 'CYS' and at.get_name() == 'SG'
    ]

    if not sg_atoms:
        print("No Cys SG atoms found in structure.")
        return

    nbsearch = NeighborSearch(sg_atoms)

    ssbonds = []
    for at1, at2 in nbsearch.search_all(args.cutoff):
        res1 = at1.get_parent()
        res2 = at2.get_parent()

        # A disulphide bond only makes sense between different Cys residues
        if res1 is res2:
            continue

        dist = at1 - at2
        ssbonds.append((res1, at1, res2, at2, dist))

    # Sort by residue number / chain of the first atom involved
    ssbonds.sort(key=lambda h: (h[0].get_parent().id, h[0].id[1], h[2].get_parent().id, h[2].id[1]))

    print(f"Possible disulphide bonds (Cys SG-SG distance < {args.cutoff} A)")
    print("-" * 70)
    for res1, at1, res2, at2, dist in ssbonds:
        print(f"{residue_id(res1):>12}.{at1.get_name():<4} -- "
              f"{residue_id(res2):<12}.{at2.get_name():<4}  {dist:6.2f} A")

    print(f"\nTotal possible disulphide bonds found: {len(ssbonds)}")


if __name__ == '__main__':
    main()