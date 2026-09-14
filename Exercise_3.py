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

print("Possible hydrogen bonds (polar atom-atom distance <", args.cutoff, "A)")
print()

ncontact = 0

for at1, at2 in nbsearch.search_all(args.cutoff):
    res1 = at1.get_parent()
    res2 = at2.get_parent()
    if res1 == res2:
        continue
    ncontact += 1
    chain1 = res1.get_parent().id
    chain2 = res2.get_parent().id
    dist = at1 - at2

    print(ncontact, res1.get_resname(), chain1 + str(res1.id[1]) + "." + at1.get_name(),
          "-", res2.get_resname(), chain2 + str(res2.id[1]) + "." + at2.get_name(),
          ": %.2f A" % dist)

print()
print("Total possible hydrogen bonds found:", ncontact)