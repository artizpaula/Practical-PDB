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
    print(f"No residue with number {args.res_num} found"
          f"{' in chain ' + args.chain if args.chain else ''}.")

for res in hits:
    chain_id = res.get_parent().id
    print(f"\nResidue: {res.get_resname()} {chain_id}{res.id[1]}")
    print("-" * 50)
    print(f"{'Atom':<6}{'X':>10}{'Y':>10}{'Z':>10}")
    atom_list = sorted(res.get_atoms(), key=lambda a: a.get_serial_number())
    for atom in atom_list:
        x, y, z = atom.get_coord()
        print(f"{atom.get_name():<6}{x:>10.3f}{y:>10.3f}{z:>10.3f}")
