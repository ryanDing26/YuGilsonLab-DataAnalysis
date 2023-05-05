#! /usr/bin/env python
# Vina Docking sample

'''
Dictionary of atomic weights of different elements
'''
ATOMIC_WEIGHTS = {'H':1.008, 'HE':4.002602, 'LI':6.94, 'BE':9.012182,
       'B':10.81, 'C':12.011, 'N':14.007, 'O':15.999, 'F':18.9984032,
       'NE':20.1797, 'NA':22.98976928, 'MG':24.305, 'AL':26.9815386,
       'SI':28.085, 'P':30.973762, 'S':32.06, 'CL':35.45, 'AR':39.948,
       'K':39.0983, 'CA':40.078, 'SC':44.955912, 'TI':47.867, 'V':50.9415,
       'CR':51.9961, 'MN':54.938045, 'FE':55.845, 'CO':58.933195,
       'NI':58.6934, 'CU':63.546, 'ZN':65.38, 'GA':69.723, 'GE':72.630,
       'AS':74.92160, 'SE':78.96, 'BR':79.904, 'RB':85.4678, 'SR':87.62,
       'Y':88.90585, 'ZR':91.224, 'NB':92.90638, 'MO':95.96, 'TC':98,
       'RU':101.07, 'RH':102.90550, 'PD':106.42, 'AG':107.8682, 'CD':112.411,
       'IN':114.818, 'SN':118.710, 'SB':121.760, 'TE':127.60, 'I':126.90447,
       'XE':131.293, 'CS':132.9054519, 'BA':137.327, 'LA':138.90547,
       'CE':140.116, 'PR':140.90765, 'ND':144.242, 'PM':145, 'SM':150.36,
       'EU':151.964, 'GD':157.25, 'TB':158.92535, 'DY':162.500, 'HO':164.93032,
       'ER':167.259, 'TM':168.93421, 'YB':173.054, 'LU':174.9668, 'HF':178.49,
       'TA':180.94788, 'W':183.84, 'RE':186.207, 'OS':190.23, 'IR':192.217,
       'PT':195.084, 'AU':196.966569, 'HG':200.592, 'TL':204.38, 'PB':207.2,
       'BI':208.98040, 'PO':209, 'AT':210, 'RN':222, 'FR':223, 'RA':226,
       'AC':227, 'TH':232.03806, 'PA':231.03588, 'U':238.02891, 'NP':237,
       'PU':244, 'AM':243, 'CM':247, 'BK':247, 'CF':251, 'ES':252, 'FM':257,
       'MD':258, 'NO':259, 'LR':262, 'RF':267, 'DB':268, 'SG':269, 'BH':270,
       'HS':269, 'MT':278, 'DS':281, 'RG':281, 'CN':285, 'UUT':286, 'FL':289,
       'UUP':288, 'LV':293, 'UUS':294}

import os
import subprocess
from vina import Vina
import sys

os.system("start vina.exe")
v = Vina(sf_name='vina', cpu = 2, seed = 12345)

def dock_vina(struct, receptor_file, ligand_file, ):
    subprocess.run( ["bash", "pdtqtconverter.bash", sys.argv[1]])
    v.set_receptor(receptor_file)
    v.set_ligand_from_file(ligand_file)
    v.compute_vina_maps(center=center_of_mass("/Users/yashravipati/Downloads/PDBBind_processed/%s/%s_protein_processed.pdb", struct, struct), box_size=[20, 20, 20])
    v.dock(exhaustiveness=32, n_poses=1, max_evals=1000000000)
    v.write_poses(ligand_file[:len("/Users/yashravipati/Downloads/PDBBind_processed/") + 9] + "_docked.pdbqt", n_poses=1, overwrite=True)

# Code pulled from https://github.com/shamsaraj/CrossDocker/blob/master/source/center.py
def center_of_mass(pdbfile, include='ATOM,HETATM'):
    """
    Calculates center of mass of a protein and/or ligand structure.
    Returns:
        center (list): List of float coordinates [x,y,z] that represent the
        center of mass (precision 3).
    """
    center = [None, None, None]
    include = tuple(include.split(','))
    
    # Issue: not sure if this works with pdbqt or just pdb files to read in coordinates;
    # If pdb, then after this we can run the bash script to convert a file to pdbqt
    with open(pdbfile, 'r') as pdb:

        # extract coordinates [ [x1,y1,z1], [x2,y2,z2], ... ]
        coordinates = []
        masses = []    
        for line in pdb:
            if line.startswith(include):
                coordinates.append([float(line[30:38]),    # x_coord
                                    float(line[38:46]),    # y_coord
                                    float(line[46:54])     # z_coord
                                   ])
                element_name = line[76:].strip()
                if element_name not in ATOMIC_WEIGHTS:
                    element_name = line.split()[2].strip()[0]
                masses.append(ATOMIC_WEIGHTS[element_name])

        assert len(coordinates) == len(masses)

        # calculate relative weight of every atomic mass
        total_mass = sum(masses)
        weights = [float(atom_mass/total_mass) for atom_mass in masses]

        # calculate center of mass
        center = [sum([coordinates[i][j] * weights[i]
              for i in range(len(weights))]) for j in range(3)]
        center_rounded = [round(center[i], 3) for i in range(3)]
        return center_rounded

def main():
    # Idea: Create a GUI that would take in a certain file of ligands instead of having these two variables preset by a user
    root_dir = "/Users/yashravipati/Downloads/PDBBind_processed"
    # Currently just takes in a file name that the user specifies
    struct = ""
    # Type exit to stop docking ligands
    while (struct != "exit"):
        struct = input("Enter the ligand here (type exit to stop): ")
        dock_vina("/Users/yashravipati/Downloads/PDBBind_processed/%s/%s_protein_processed.pdb.pdbqt",
              "/Users/yashravipati/Downloads/PDBBind_processed/%s/%s_ligand.mol2.pdbqt", struct, struct, struct, struct)
        
if __name__ == "__main__": main()