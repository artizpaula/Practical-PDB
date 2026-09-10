# Practical 1 - PDB

# Exercise 1

#!/usr/bin/env python

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

    # Crear el parser de PDB
    parser_pdb = PDBParser(PERMISSIVE=1)

    # Cargar la estructura desde el archivo PDB
    st = parser_pdb.get_structure('structure', args.pdb_file)

    select = []

    # Seleccionar únicamente los átomos CA
    print("Átomos CA seleccionados:")

    for at in st.get_atoms():
        if at.id == 'CA':
            select.append(at)

            print(
                f"ATOM: {at.get_parent().get_resname()}, "
                f"{at.get_parent().id[1]}, "
                f"{at.id}"
            )

    # Preparar la búsqueda de vecinos
    nbsearch = NeighborSearch(select)

    print("\nNBSEARCH:")
    print(f"Distancia máxima: {args.distance} Å\n")

    # Buscar contactos
    ncontact = 1

    for at1, at2 in nbsearch.search_all(args.distance):
        print(f"Contact: {ncontact}")

        print(
            f"at1: {at1}, "
            f"{at1.get_serial_number()}, "
            f"{at1.get_parent().get_resname()}"
        )

        print(
            f"at2: {at2}, "
            f"{at2.get_serial_number()}, "
            f"{at2.get_parent().get_resname()}"
        )

        print()

        ncontact += 1

    print(f"Total de contactos encontrados: {ncontact - 1}")


if __name__ == '__main__':
    main()