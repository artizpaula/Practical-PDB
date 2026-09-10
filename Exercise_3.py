# Exercise 3

#!/usr/bin/env python

import argparse
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser


def main():

    # Configuración de los argumentos
    parser = argparse.ArgumentParser(
        prog='ProgName',
        description='Search for possible hydrogen bonds between polar atoms'
    )

    parser.add_argument(
        '--PDB',
        dest='pdb_file',
        help='PDB file to analyze',
        required=True
    )

    parser.add_argument(
        '--distance',
        dest='distance',
        type=float,
        default=3.5,
        help='Maximum distance for a possible hydrogen bond (default: 3.5 Å)'
    )

    args = parser.parse_args()

    # Crear el parser de PDB
    parser_pdb = PDBParser(PERMISSIVE=1)

    # Cargar la estructura
    st = parser_pdb.get_structure('structure', args.pdb_file)

    selected = []

    # Seleccionar átomos polares: O, N y S
    polar_atoms = ['O', 'N', 'S']

    print("Polar atoms selected:")

    for atom in st.get_atoms():

        # Obtener el nombre del átomo
        atom_name = atom.get_name().strip()

        if atom_name.startswith(tuple(polar_atoms)):
            selected.append(atom)

            residue = atom.get_parent()

            print(
                f"ATOM: {residue.get_resname()}, "
                f"Residue: {residue.id[1]}, "
                f"Atom: {atom_name}"
            )

    # Preparar la búsqueda de vecinos
    nbsearch = NeighborSearch(selected)

    print("\nPossible Hydrogen Bonds:")
    print(f"Distance criterion: < {args.distance} Å\n")

    ncontact = 1

    # Buscar pares de átomos próximos
    for atom1, atom2 in nbsearch.search_all(args.distance):

        # Calcular la distancia entre los dos átomos
        distance = atom1 - atom2

        # Obtener información de los residuos
        residue1 = atom1.get_parent()
        residue2 = atom2.get_parent()

        print(f"Contact: {ncontact}")

        print(
            f"Atom 1: {atom1.get_name()}, "
            f"Residue: {residue1.get_resname()} "
            f"{residue1.id[1]}, "
            f"Chain: {residue1.get_parent().id}"
        )

        print(
            f"Atom 2: {atom2.get_name()}, "
            f"Residue: {residue2.get_resname()} "
            f"{residue2.id[1]}, "
            f"Chain: {residue2.get_parent().id}"
        )

        print(f"Distance: {distance:.2f} Å")
        print()

        ncontact += 1

    print(f"Total possible hydrogen bonds: {ncontact - 1}")


if __name__ == '__main__':
    main()