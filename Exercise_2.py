# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 2

import argparse
from Bio.PDB.PDBParser import PDBParser

# argparse, following ex_cmd_line.py
parser = argparse.ArgumentParser(
    prog='Exercise_2',
    description='List all atoms (name and coordinates) for a given residue number')

parser.add_argument(
    '--chain',
    dest='chain',
    default=None,
    help='Chain id (optional). If not given, all chains are searched')

parser.add_argument(
    'input_pdb',
    help='Input PDB file')

parser.add_argument(
    'res_num',
    type=int,
    help='Residue number')

args = parser.parse_args()

pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)

# load structure from PDB file (as in ex_list_res.py)
st = pdb_parser.get_structure('structure', args.input_pdb)

model0 = st[0]

hits = []
for chain in model0:
    if args.chain is not None and chain.id != args.chain:
        continue
    for res in chain:
        if res.id[1] == args.res_num:
            hits.append(res)

if not hits:
    print("No residue with number", args.res_num, "found.")

for res in hits:
    chain_id = res.get_parent().id
    print()
    print("Residue:", res.get_resname(), chain_id + str(res.id[1]))
    atom_list = sorted(res.get_atoms(), key=lambda a: a.get_serial_number())
    for atom in atom_list:
        x, y, z = atom.get_coord()
        print("Atom:", atom.get_name(), " X:", "%.3f" % x, " Y:", "%.3f" % y, " Z:", "%.3f" % z)