#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# My first example with AutoDock Vina in python
#
import os
import subprocess
from vina import Vina


v = Vina(sf_name='vina', cpu = 6, seed = 12345)


def dock_vina(receptor_file, ligand_file):
    v.set_receptor(receptor_file)
    v.set_ligand_from_file(ligand_file)
    v.compute_vina_maps(center=[0, 0, 0], box_size=[20, 20, 20])
    v.dock(exhaustiveness=32, n_poses=1, max_evals=1000000000)
    v.write_poses(ligand_file[:len("/Users/yashravipati/Downloads/PDBBind_processed/") + 9] + "_docked.pdbqt", n_poses=1, overwrite=True)

root_dir = "/Users/yashravipati/Downloads/PDBBind_processed"

struct = "6s9c"

dock_vina("/Users/yashravipati/Downloads/PDBBind_processed/" + struct + "/" + struct + "_protein_processed.pdb.pdbqt", "/Users/yashravipati/Downloads/PDBBind_processed/" + struct + "/" + struct + "_ligand.mol2.pdbqt")