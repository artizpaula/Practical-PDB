# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 5
"""
Usage (from terminal):
    python Exercise_5.py structure.pdb
    python Exercise_5.py structure.pdb --cutoff 2.2
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
        prog='Exercise_5',
        description='List possible peptide bonds (C-N contacts) between residues')
 
    argp.add_argument(
        '--cutoff',
        dest='cutoff',
        type=float,
        default=2.0,
        help='Distance criterium (dist < 2.0 A)')
 
    argp.add_argument(
        'input_pdb',
        help='PDB file')
 
    opts = argp.parse_args()
 
    struct_name = os.path.splitext(os.path.basename(opts.input_pdb))[0]
 
    reader = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = reader.get_structure(struct_name, opts.input_pdb)
 
    backbone_atoms = [a for a in structure.get_atoms() if a.get_name() in ('C', 'N')]
 
    searcher = NeighborSearch(backbone_atoms)
 
    links = []
    for a1, a2 in searcher.search_all(opts.cutoff):
        tags = {a1.get_name(), a2.get_name()}
        if tags != {'C', 'N'}:
            continue
 
        if a1.get_name() != 'C':
            a1, a2 = a2, a1
 
        r1 = a1.get_parent()
        r2 = a2.get_parent()
 
        if r1 is r2:
            continue
 
        d = a1 - a2
        links.append((r1, a1, r2, a2, d))
 
    links.sort(key=lambda h: (h[0].get_parent().id, h[0].id[1], h[2].get_parent().id, h[2].id[1]))
 
    print(f"Possible peptide bonds (C-N distance < {opts.cutoff} A)")
    print("-" * 70)
    for r1, a1, r2, a2, d in links:
        print(f"{get_label(r1):>12}.{a1.get_name():<4} -- "
              f"{get_label(r2):<12}.{a2.get_name():<4}  {d:6.2f} A")
 
    print(f"\nTotal possible peptide bonds found: {len(links)}")
 
 
if __name__ == '__main__':
    main()
 