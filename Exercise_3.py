# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 3

import argparse
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

polar_atoms = ('O', 'N', 'S')

# argparse, following ex_cmd_line.py
parser = argparse.ArgumentParser(
    prog='Exercise_3',
    description='List possible hydrogen bonds between polar atoms')

parser.add_argument(
    '--cutoff',
    dest='cutoff',
    type=float,
    default=3.5,
    help='Distance criterium in Angstrom (default: 3.5)')

parser.add_argument(
    'input_pdb',
    help='Input PDB file')

args = parser.parse_args()

pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)

# load structure from PDB file
st = pdb_parser.get_structure('structure', args.input_pdb)

select = []

# only polar atoms (O, N, S)
for at in st.get_atoms():
    if at.element in polar_atoms:
        select.append(at)

nbsearch = NeighborSearch(select)

print(f"Possible hydrogen bonds (polar atom-atom distance < {args.cutoff} A)")

ncontact = 1

for at1, at2 in nbsearch.search_all(args.cutoff):
    res1 = at1.get_parent()
    res2 = at2.get_parent()
    if res1 == res2:
        continue

    print(f"Contact {ncontact}:")
    print(f"  {res1.get_resname()} {res1.get_parent().id}{res1.id[1]}.{at1.get_name()}"
          f" -- {res2.get_resname()} {res2.get_parent().id}{res2.id[1]}.{at2.get_name()}")
    print(f"  Distance: {at1 - at2:.2f} A")
    print()
    ncontact += 1
print(f"Total possible hydrogen bonds found: {ncontact - 1}")