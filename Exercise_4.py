# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 4

import argparse
 
from Bio.PDB.PDBParser import PDBParser
 
AA3_TO_AA1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}
AA1_TO_AA3 = {one: three for three, one in AA3_TO_AA1.items()}
 
def to_three_letter(code):
    code = code.strip().upper()
    if len(code) == 1:
        if code not in AA1_TO_AA3:
            raise SystemExit(f"Unknown one-letter residue code: {code}")
        return AA1_TO_AA3[code]
    elif len(code) == 3:
        if code not in AA3_TO_AA1:
            raise SystemExit(f"Unknown three-letter residue code: {code}")
        return code
    else:
        raise SystemExit(f"Residue type must be a 1- or 3-letter code, got: {code}")
 
 
# argparse, following ex_cmd_line.py
parser = argparse.ArgumentParser(
    prog='Exercise_4',
    description='List all CA atoms (with coordinates) of a given residue type')
 
parser.add_argument(
    'input_pdb',
    help='Input PDB file')
 
parser.add_argument(
    'res_type',
    help='Residue type, either one-letter (e.g. R) or three-letter (e.g. ARG) code')
 
args = parser.parse_args()
 
target_resname = to_three_letter(args.res_type)
 
pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)
 
# load structure from PDB file (as in ex_list_res.py)
st = pdb_parser.get_structure('structure', args.input_pdb)
 
matches = []
for res in st.get_residues():
    if res.get_resname() == target_resname:
        for at in res:
            if at.id == 'CA':
                matches.append(at)
 
matches.sort(key=lambda a: (a.get_parent().get_parent().id, a.get_parent().id[1]))
 
print(f"CA atoms of residue type {target_resname} ({AA3_TO_AA1[target_resname]})")
print(f"{'Chain':<6}{'ResNum':<8}{'X':>10}{'Y':>10}{'Z':>10}")
for at in matches:
    res = at.get_parent()
    chain_id = res.get_parent().id
    x, y, z = at.get_coord()
    print(f"{chain_id:<6}{res.id[1]:<8}{x:>10.3f}{y:>10.3f}{z:>10.3f}")
 
print(f"\nTotal {target_resname} residues found: {len(matches)}")
