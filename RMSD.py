import requests
import os
import subprocess
from Bio.PDB import *
import numpy as np
import rmsd

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Superimposer import Superimposer

def calculate_rmsd(native_pdb, rcsb_pdb):
    parser = PDBParser(QUIET=True)
    native_structure = parser.get_structure("native", native_pdb)
    rcsb_structure = parser.get_structure("rcsb", rcsb_pdb)

    super_imposer = Superimposer()
    super_imposer.set_atoms(native_structure.get_atoms(), rcsb_structure.get_atoms())
    super_imposer.apply(rcsb_structure.get_atoms())

    rmsd = super_imposer.rms
    return rmsd


def rmsd(structure1, structure2):
    super_imposer = Superimposer()
    super_imposer.set_atoms(structure1.get_atoms(), structure2.get_atoms())
    return super_imposer.rms


def get_structure(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        return None

root_dir = "/Users/yashravipati/Downloads/DockedFiles"

file = open("/Users/yashravipati/Documents/rmsd_values.txt", "w")

print(calculate_rmsd())
x = 0
for dirpath, dirnames, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.endswith(".pdb"):
            pdb_id = filename[:4]
            structurePath = root_dir + "/" + filename
            structureRCSBPath = get_structure(pdb_id)

            print (structureRCSBPath)

            rmsdval = calculate_rmsd(structurePath, structureRCSBPath)
            file.write("RMSD value: " + str(rmsdval) + "\n")
            print (rmsdval)

file.close()