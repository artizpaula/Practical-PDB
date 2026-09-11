# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python
#
""" Exercise 3
Determine all possible hydrogen bonds (Polar atoms at less than 3.5 A).

Parameters: PDB file name. Optional: cut-off distance (defaults to 3.5)

Usage (from terminal):
    python Exercise_3.py structure.pdb
    python Exercise_3.py structure.pdb --cutoff 3.2
"""

import argparse
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

POLAR_ELEMENTS = ('O', 'N', 'S')


def residue_id(res):
    chain_id = res.get_parent().id
    resnum = res.id[1]
    icode = res.id[2].strip()
    return f"{res.get_resname()} {chain_id}{resnum}{icode}"


def main():
    parser = argparse.ArgumentParser(
        prog='Exercise_3',
        description='List possible hydrogen bonds between polar atoms (O, N, S) closer than a cut-off distance'
    )

    parser.add_argument(
        '--cutoff',
        dest='cutoff',
        type=float,
        default=3.5,
        help='Distance cut-off in Angstroms for a possible hydrogen bond (default: 3.5)'
    )

    parser.add_argument(
        'pdb_file',
        help='Input PDB file'
    )

    args = parser.parse_args()

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]

    pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)
    st = pdb_parser.get_structure(pdb_id, args.pdb_file)

    # Select only polar atoms (O, N, S)
    polar_atoms = [at for at in st.get_atoms() if at.element in POLAR_ELEMENTS]

    nbsearch = NeighborSearch(polar_atoms)

    hbonds = []
    for at1, at2 in nbsearch.search_all(args.cutoff):
        res1 = at1.get_parent()
        res2 = at2.get_parent()
        # Hydrogen bonds only make sense between atoms of different residues
        if res1 is res2:
            continue
        dist = at1 - at2
        hbonds.append((res1, at1, res2, at2, dist))

    # Sort by residue number / chain of the first atom involved
    hbonds.sort(key=lambda h: (h[0].get_parent().id, h[0].id[1], h[2].get_parent().id, h[2].id[1]))

    print(f"Possible hydrogen bonds (polar atom-atom distance < {args.cutoff} A)")
    print("-" * 70)
    for res1, at1, res2, at2, dist in hbonds:
        print(f"{residue_id(res1):>12}.{at1.get_name():<4} -- "
              f"{residue_id(res2):<12}.{at2.get_name():<4}  {dist:6.2f} A")

    print(f"\nTotal possible hydrogen bonds found: {len(hbonds)}")


if __name__ == '__main__':
    main()