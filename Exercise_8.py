# Practical 1 - PDB
# Alicia Mañas, Lídia Sanchez and Paula Artiz

#!/usr/bin/env python

# Exercise 8
from Bio.PDB import PDBList
import os

def ex_8(id):
    # Strip any extension that maybe has been included that affects to find the id
    pdb_id = os.path.splitext(id)[0].lower()

    new_filename = f"{pdb_id}.pdb"

    pdb_list = PDBList()
    pdb_list.retrieve_pdb_file(
        pdb_id,
        pdir=".",
        file_format="pdb"
    )

    old_filename = f"./pdb{pdb_id}.ent" #retrieve_pdb_file always save the pdb file with the .ent extension

    os.rename(old_filename, new_filename)

    return new_filename
    
