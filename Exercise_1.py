# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 1

import argparse
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

# argparse, following ex_cmd_line.py
parser = argparse.ArgumentParser(
    prog='Exercise_1',
    description='List pairs of residues whose CA atoms are closer than a given distance')

parser.add_argument(
    'input_pdb',
    help='Input PDB file')

parser.add_argument(
    'distance',
    type=float,
    help='Distance threshold (Angstrom) to consider two residues in contact')

args = parser.parse_args()

pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)

# load structure from PDB file (as in ex_distances.py)
st = pdb_parser.get_structure('structure', args.input_pdb)

select = []

# Select only CA atoms
for at in st.get_atoms():
    if at.id == 'CA':
        select.append(at)

# Preparing search
nbsearch = NeighborSearch(select)

print(f"Pairs of residues with CA atoms closer than {args.distance} A")

# Searching for contacts under the given distance
ncontact = 1

for at1, at2 in nbsearch.search_all(args.distance):
    res1 = at1.get_parent()
    res2 = at2.get_parent()

    if res1 == res2:
        continue

    print(f"Contact {ncontact}:")
    print(f"  Residue 1: {res1.get_resname()} {res1.get_parent().id}{res1.id[1]}"
          f"  (CA serial {at1.get_serial_number()})")
    print(f"  Residue 2: {res2.get_resname()} {res2.get_parent().id}{res2.id[1]}"
          f"  (CA serial {at2.get_serial_number()})")
    print(f"  Distance: {at1 - at2:.2f} A")
    print()
    ncontact += 1

print(f"Total pairs found: {ncontact - 1}")