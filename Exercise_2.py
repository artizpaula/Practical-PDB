# Exercise 2

""" Simple program to print ARG residues iteration over atoms """

import argparse
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

def main():
    # Configuración de los argumentos de línea de comandos
    parser = argparse.ArgumentParser(
        prog='ProgName',
        description='Programa para buscar contactos entre átomos CA de una estructura PDB'
    )

    parser.add_argument(
        '--PDB',
        dest='pdb_file',
        help='Archivo PDB que se quiere analizar',
        required=True
    )

    parser.add_argument(
        '--distance',
        dest='distance',
        type=float,
        default=20.0,
        help='Distancia máxima para considerar un contacto (por defecto: 20 Å)'
    )

    args = parser.parse_args()

parser = PDBParser()

st = parser.get_structure('1UBQ', '1ubq.pdb')

selected = []
aa = ["ARG"]

for at in st.get_atoms():
    if at.get_parent().get_resname() in aa:
        selected.append(at)

print("Coordinates:")
for atom in selected:
    print(f"{atom.get_parent().get_resname()}, {atom.get_parent().id}, {atom.get_name()}, {atom.get_coord()}")

if __name__ == '__main__':
    main()