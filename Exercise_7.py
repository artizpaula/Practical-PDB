# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 7

from Bio.PDB.PDBParser import PDBParser
import numpy as np
import argparse
parser = argparse.ArgumentParser(
    prog='Exercise_7',
    description='Print distances between all atom pairs of two given residues')

parser.add_argument(
    'input_pdb',
    help='Input PDB file')
parser.add_argument(
    'residue_1',
    help='Introduce residue 1'
)
parser.add_argument(
    'residue_2',
    help='Introduce residue 2'
)
args = parser.parse_args()

pdb_parser = PDBParser() #The way that we can read the pdb
st = pdb_parser.get_structure('protein', args.input_pdb)


res1 = st[0]["A"][args.residue_1]
#Selection of Residues 1 and 2of Chain A
res2 = st[0]["A"][args.residue_2]

print("Residue 1 is", res1.get_resname())
print("Residue 2 is", res2.get_resname())

print("\nAtom1 Atom2 dist1 dist2\n-------------------------")
for at1 in res1.get_atoms():      # Replace get_atoms with get_atom if you get an Error!
    for at2 in res2.get_atoms():
        dist = at2 - at1    # Direct procedure with (-) to compute distances
        vector = at2.coord - at1.coord  # Or using numpy coordinates
        distance = np.sqrt(np.sum(vector ** 2))
        print(at1, at2, dist, distance)

