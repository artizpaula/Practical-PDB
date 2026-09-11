# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

""" Exercise 1
Determine the list of pairs of residues whose CA atoms are closer than a
given distance.

Parameters: PDB file name, distance.

Usage (from terminal):
    python Exercise_1.py structure.pdb 5.0
    python Exercise_1.py --dist 5.0 structure.pdb
"""

import argparse
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch


def residue_id(res):
    """Return a readable identifier for a residue: RESNAME ChainId ResNum"""
    chain_id = res.get_parent().id
    resnum = res.id[1]
    icode = res.id[2].strip()
    return f"{res.get_resname()} {chain_id}{resnum}{icode}"


def main():
    parser = argparse.ArgumentParser(
        prog='Exercise_1',
        description='List pairs of residues whose CA atoms are closer than a given distance'
    )

    parser.add_argument(
        '--dist',
        dest='dist',
        type=float,
        default=5.0,
        help='Distance threshold in Angstroms (default: 5.0)'
    )

    parser.add_argument(
        'pdb_file',
        help='Input PDB file'
    )

    args = parser.parse_args()

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]

    pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)
    st = pdb_parser.get_structure(pdb_id, args.pdb_file)

    # Select only CA atoms
    ca_atoms = [at for at in st.get_atoms() if at.id == 'CA']

    nbsearch = NeighborSearch(ca_atoms)

    pairs = []
    for at1, at2 in nbsearch.search_all(args.dist):
        res1 = at1.get_parent()
        res2 = at2.get_parent()
        # Skip pairs that are the same residue (should not happen for CA-CA
        # but kept as a safety check) and avoid trivial neighbor pairs
        if res1 is res2:
            continue
        dist = at1 - at2
        pairs.append((res1, res2, dist))

    # Sort by residue number of the first residue, then the second
    pairs.sort(key=lambda p: (p[0].get_parent().id, p[0].id[1], p[1].get_parent().id, p[1].id[1]))

    print(f"Pairs of residues with CA-CA distance < {args.dist} A")
    print("-" * 60)
    for res1, res2, dist in pairs:
        print(f"{residue_id(res1):>12} -- {residue_id(res2):<12}  {dist:6.2f} A")

    print(f"\nTotal pairs found: {len(pairs)}")


if __name__ == '__main__':
    main()