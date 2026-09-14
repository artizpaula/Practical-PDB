# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 6

"""
Usage (from terminal):
    python Exercise_6.py structure.pdb
    python Exercise_6.py structure.pdb --cutoff 2.2
"""

import argparse
import os
 
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
 
 
def get_label(residue):
    seg = residue.get_parent().id
    num = residue.id[1]
    ins = residue.id[2].strip()
    return f"{residue.get_resname()} {seg}{num}{ins}"
 
 
def main():
    argp = argparse.ArgumentParser(
        prog='Exercise_6',
        description='List possible disulphide bonds (Cys SG-SG contacts)')
 
    argp.add_argument(
        '--cutoff',
        dest='cutoff',
        type=float,
        default=1.9,
        help='Distance criterium (dist < 1.9 A)')
 
    argp.add_argument(
        'input_pdb',
        help='PDB file')
 
    opts = argp.parse_args()
 
    struct_name = os.path.splitext(os.path.basename(opts.input_pdb))[0]
 
    reader = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = reader.get_structure(struct_name, opts.input_pdb)
 
    sg_list = [a for a in structure.get_atoms()
        if a.get_parent().get_resname() == 'CYS' and a.get_name() == 'SG']
 
    if not sg_list:
        print("No Cys SG atoms found in structure.")
        return
 
    searcher = NeighborSearch(sg_list)
 
    bridges = []
    for a1, a2 in searcher.search_all(opts.cutoff):
        r1 = a1.get_parent()
        r2 = a2.get_parent()
 
        if r1 is r2:
            continue
 
        d = a1 - a2
        bridges.append((r1, a1, r2, a2, d))
 
    bridges.sort(key=lambda h: (h[0].get_parent().id, h[0].id[1], h[2].get_parent().id, h[2].id[1]))
 
    print(f"Possible disulphide bonds (Cys SG-SG distance < {opts.cutoff} A)")
    print("-" * 70)
    for r1, a1, r2, a2, d in bridges:
        print(f"{get_label(r1):>12}.{a1.get_name():<4} -- "
              f"{get_label(r2):<12}.{a2.get_name():<4}  {d:6.2f} A")
 
    print(f"\nTotal possible disulphide bonds found: {len(bridges)}")
 
 
if __name__ == '__main__':
    main()