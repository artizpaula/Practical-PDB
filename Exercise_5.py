 
import argparse
import os
 
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
 
 
def residue_id(res):  # residue as ASP A32
    return f"{res.get_resname()} {res.get_parent().id}{str(res.id[1])}"
 
 
def atom_id(atom):  # atom as ASP A32.N, re-uses residue_id
    return f"{residue_id(atom.get_parent())}.{atom.id}"
 
 
parser = argparse.ArgumentParser(
    prog='Exercise_5',
    description='Generate a list of backbone connectivity (peptide bonds, C-N contacts)'
)
 
parser.add_argument(
    '--cutoff',
    dest='cutoff',
    type=float,
    default=2.5,
    help='Cut-off distance for peptide bonds (default: 2.5 A)'
)
 
parser.add_argument(
    'pdb_file',
    help='PDB file'  # 1UBQ.pdb or 4HHB.pdb
)
 
args = parser.parse_args()
 
print("\nSettings\n--------")
for k, v in vars(args).items():
    print('{:10}:'.format(k), v)
print()
 
CUTOFF = args.cutoff
pdb_file = args.pdb_file
pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
 
 
pdb_parser = PDBParser(PERMISSIVE=1)
 
st = pdb_parser.get_structure(pdb_id, pdb_file)
 
# Select only backbone C and N atoms (identified by atom id, not element)
 
select = []
 
for at in st.get_atoms():
    if at.id in ('C', 'N'):
        select.append(at)
 
print(f"{len(select)} backbone C/N atoms selected")
 
nbsearch = NeighborSearch(select)
 
print("\nNBSEARCH:")

 
pbonds = []
ncontact = 1
 
for at1, at2 in nbsearch.search_all(CUTOFF):
 
    # Only interested in a C paired with a N (not C-C or N-N)
    ids = {at1.id, at2.id}
    if ids != {'C', 'N'}:
        continue
 
    # Make sure at1 is always the C and at2 is always the N
    if at1.id != 'C':
        at1, at2 = at2, at1
 
    res1 = at1.get_parent()
    res2 = at2.get_parent()
 
    # A peptide bond only makes sense between different residues
    if res1 is res2:
        continue
 
    dist = at1 - at2  # Direct procedure with (-) to compute distances
 
    pbonds.append((res1, at1, res2, at2, dist))
 
    print(f"Contact: {ncontact}")
    print(f"at1: {atom_id(at1)}, {at1.get_serial_number()}")
    print(f"at2: {atom_id(at2)}, {at2.get_serial_number()}")
    print(f"dist: {dist:.2f} A")
    print()
 
    ncontact += 1

print(f"Total peptide bonds found: {len(pbonds)}")