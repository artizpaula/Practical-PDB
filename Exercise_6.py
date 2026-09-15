
import argparse
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch


def residue_id(res):  # residue as ASP A32
    return f"{res.get_resname()} {res.get_parent().id}{str(res.id[1])}"


def atom_id(atom):  # atom as ASP A32.N, re-uses residue_id
    return f"{residue_id(atom.get_parent())}.{atom.id}"


parser = argparse.ArgumentParser(
    prog='Exercise_6',
    description='Generate a list of disulphide bond connectivity (Cys SG-SG contacts)'
)

parser.add_argument(
    '--cutoff',
    dest='cutoff',
    type=float,
    default=2.5,
    help='Distance criterium (default: 2.5 A)'
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


pdb_parser = PDBParser(PERMISSIVE=1, QUIET = True)

st = pdb_parser.get_structure(pdb_id, pdb_file)

# Select only the side-chain sulfur (SG) atoms of Cys residues

select = []

for at in st.get_atoms():
    if at.get_parent().get_resname() == 'CYS' and at.id == 'SG':
        select.append(at)

print(f"{len(select)} Cys SG atoms selected")

if not select:
    print("\nNo Cys SG atoms found in structure.")
    raise SystemExit(0)

nbsearch = NeighborSearch(select)

print("\nNBSEARCH:")

ssbonds = []
ncontact = 1

for at1, at2 in nbsearch.search_all(CUTOFF):

    res1 = at1.get_parent()
    res2 = at2.get_parent()

    # A disulphide bond only makes sense between different Cys residues
    if res1 is res2:
        continue

    dist = at1 - at2  

    ssbonds.append((res1, at1, res2, at2, dist))

    print(f"Contact: {ncontact}")
    print(f"at1: {atom_id(at1)}, {at1.get_serial_number()}")
    print(f"at2: {atom_id(at2)}, {at2.get_serial_number()}")
    print(f"dist: {dist:.2f} A")
    print()

    ncontact += 1

print(f"Total disulphide bonds found: {len(ssbonds)}")
