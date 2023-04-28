import subprocess
# tqdm is literally just a progress bar
from tqdm import tqdm
import os
# RDKit packages
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem import MolFromSmiles, QED
from rdkit.Chem import AllChem
# Synthetic Accessibility Scorer
from sascorer import calculateScore
import numpy as np
import time
import math
import json
import torch
import csv
import pickle

delta_g_to_kd = lambda x: math.exp(x / (0.00198720425864083 * 298.15))
kd_to_delta_g = lambda x: 0.00198720425864083 * 298.15 * math.log(x)
abfe_devices = [5, 6, 7]
proc_max = 16
# surrogate_model = pickle.load(open('docking_surrogate.pkl', 'rb'))
 
def smiles_to_sa(smiles):
    vals = []
    for smile in tqdm(smiles):
        vals.append(calculateScore(MolFromSmiles(smile)))
    return vals

def smiles_to_qed(smiles):
    vals = []
    for smile in tqdm(smiles):
        vals.append(QED.qed(MolFromSmiles(smile)))
    return vals

def smiles_to_logp(smiles):
    vals = []
    for smile in tqdm(smiles):
        vals.append(MolLogP(MolFromSmiles(smile)))
    return vals

def smiles_to_morgan(smiles):
    out = []
    for smile in tqdm(smiles):
        out.append(np.array(AllChem.GetMorganFingerprintAsBitVect(MolFromSmiles(smile), 3, nBits=2048)))
    return np.array(out)

def smiles_to_affinity(smiles, autodock='~/AutoDock-GPU/bin/autodock_gpu_128wi', protein_file='/data/peter/pdbs/b65d866690c8eaf8.maps.fld', num_devices=torch.cuda.device_count(), starting_device=0):
    if not os.path.exists('ligands'):
        os.mkdir('ligands')
    if not os.path.exists('outs'):
        os.mkdir('outs')
    subprocess.run('rm core.*', shell=True, stderr=subprocess.DEVNULL)
    subprocess.run('rm outs/*.xml', shell=True, stderr=subprocess.DEVNULL)
    subprocess.run('rm outs/*.dlg', shell=True, stderr=subprocess.DEVNULL)
    subprocess.run('rm -rf ligands/*', shell=True, stderr=subprocess.DEVNULL)
    for device in range(starting_device, starting_device + num_devices):
        os.mkdir(f'ligands/{device}')
    device = starting_device
    for i, smile in enumerate(tqdm(smiles, desc='preparing ligands')):
        subprocess.Popen(f'obabel -:"{smile}" -O ligands/{device}/ligand{i}HASH{hash(smile)}.pdbqt -p 7.4 --partialcharge gasteiger --gen3d', shell=True, stderr=subprocess.DEVNULL)
        device += 1
        if device == starting_device + num_devices:
            device = starting_device
    while True:
        total = 0
        for device in range(starting_device, starting_device + num_devices):
            total += len(os.listdir(f'ligands/{device}'))
        if total == len(smiles):
            break
    time.sleep(1)
    print('running autodock..')
    subprocess.run('rm outs/*.xml', shell=True, stderr=subprocess.DEVNULL)
    subprocess.run('rm outs/*.dlg', shell=True, stderr=subprocess.DEVNULL)
    if len(smiles) == 1:
        subprocess.run(f'{autodock} -M {protein_file} -L ligands/0/ligand0.pdbqt -N outs/ligand0', shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        ps = []
        for device in range(starting_device, starting_device + num_devices):
            ps.append(subprocess.Popen(f'{autodock} -M {protein_file} -B ligands/{device}/ligand*.pdbqt -N ../../outs/ -D {device + 1}', shell=True, stdout=subprocess.DEVNULL))
        stop = False
        while not stop: 
            stop = True
            for p in ps:
                if p.poll() is None:
                    time.sleep(1)
                    stop = False
    affins = [0 for _ in range(len(smiles))]
    for file in os.listdir('outs'):
        if file.endswith('.dlg'):
            content = open(f'outs/{file}').read()
            if '0.000   0.000   0.000  0.00  0.00' not in content:
                try:
                    affins[int(file.split('ligand')[1].split('HASH')[0])] = float([line for line in content.split('\n') if 'RANKING' in line][0].split()[3])
                except:
                    pass
    return [min(affin, 0) for affin in affins]

def autodock(smiles):
    affins = np.array(smiles_to_affinity(smiles * 10))
    affin_mins = affins.reshape((-1, len(smiles))).min(0)
    smile_to_data = {smiles[i]: {'total_energy': affin_mins[i]} for i in range(len(smiles))}
    for smile in smile_to_data:
        for f_name in [f for f in os.listdir('outs') if ('.dlg' in f and str(hash(smile)) in f)]:
            f = open(f'outs/{f_name}', 'r').read()
            if f"Estimated Free Energy of Binding    =  {smile_to_data[smile]['total_energy']}" in f or f"Estimated Free Energy of Binding    = {smile_to_data[smile]['total_energy']}" in f:
                if f"Estimated Free Energy of Binding    =  {smile_to_data[smile]['total_energy']}" in f:
                    pdb = f.split(f"Estimated Free Energy of Binding    =  {smile_to_data[smile]['total_energy']}")[1]
                else:
                    pdb = f.split(f"Estimated Free Energy of Binding    = {smile_to_data[smile]['total_energy']}")[1]
                smile_to_data[smile]['intermolecular_energy'] = float(f.split('(1) Final Intermolecular Energy')[1].split()[1].strip())
                smile_to_data[smile]['internal_energy'] = float(f.split('(2) Final Total Internal Energy')[1].split()[1].strip())
                smile_to_data[smile]['torsional_energy'] = float(f.split('(3) Torsional Free Energy')[1].split()[1].strip())
                smile_to_data[smile]['unbound_energy'] = float(f.split('(4) Unbound System\'s Energy')[1].split()[1].strip())
                pdb = pdb.split('DOCKED: REMARK                         _______ _______ _______ _____ _____    ______ ____')[1].split('DOCKED: ENDMDL')[0].strip()
                pdb = pdb.replace('DOCKED: ', '')
                atoms = []
                for line in pdb.split('\n'):
                    if line.startswith('ATOM'):
                        _, _, type, _, _, x, y, z, vdw, _, _, _ = line.split()
                        atoms.append((type, float(vdw), float(x), float(y), float(z)))
                smile_to_data[smile]['atom_coords'] = atoms
                smile_to_data[smile]['number_of_atoms'] = int(f.split('Number of atoms:')[1].split()[0].strip())
                smile_to_data[smile]['number_of_rotatable_bonds'] = int(f.split('Number of rotatable bonds:')[1].split()[0].strip())
                pdb += '\nENDMDL'
                new_pdb = ''
                for line in pdb.split('\n'):
                    new_pdb += line[:66] + '\n'
                pdb = new_pdb
                open('autodock_pose.pdb', 'w').write(pdb)
                # runs bash commands inside the script
                subprocess.run('pymol -cq pymol_script.py', shell=True)
                subprocess.call(f'obabel autodock_pose.pdb -O autodock_pose.pdbqt -p 7.4', shell=True, stdout=subprocess.DEVNULL)
                # Need to change this line to fit our own directories
                subprocess.call('python /home/peter/limo/binana/python/run_binana.py -receptor /data/peter/pdbs/b65d866690c8eaf8.pdbqt -ligand autodock_pose.pdbqt -output_json binana_out.json', shell=True, stdout=subprocess.DEVNULL)
                out = json.load(open('binana_out.json', 'r'))
                smile_to_data[smile]['hydrogen_bonds'] = len(out['hydrogenBonds'])
                smile_to_data[smile]['pi_pi_stacking_interactions'] = len(out['piPiStackingInteractions'])
                smile_to_data[smile]['salt_bridges'] = len(out['saltBridges'])
                smile_to_data[smile]['t_stacking_interactions'] = len(out['tStackingInteractions'])
                break
    return smile_to_data

def abfe(smiles):
    pass

def load_bindingdb_data(file):
    out = []
    for row in csv.reader(open(file, 'r'), delimiter='	'):
        if row[9] and '<' not in row[9] and '>' not in row[9] and row[9] != 'IC50 (nM)':
            out.append([row[1], kd_to_delta_g(float(row[9]) / 1e9)])
    return out

# def dock6_flex(smiles):
#     smile_to_data = {}
#     for start_index in range(0, len(smiles), proc_max):
#         subprocess.run('rm -rf dock6/runs', shell=True, stderr=subprocess.DEVNULL)
#         os.mkdir('dock6/runs')
#         ps = []
#         for i, smile in enumerate(smiles[start_index:start_index + proc_max]):
#             dir = f'dock6/runs/{i}'
#             os.mkdir(dir)
#             ps.append(subprocess.Popen(f'obabel -:"{smile}" -O {dir}/ligand.mol2 -h --partialcharge mmff94 --gen3d', shell=True))
#         for p in ps:
#             p.wait()
#         open('run_dock6_flex_modified', 'w').write(open('run_dock6_flex', 'r').read().replace('_MAX_FILE_', str(len(smiles[start_index:start_index + proc_max]) - 1)))
#         subprocess.run(f'./run_dock6_flex_modified')
#         for i, smile in enumerate(smiles[start_index:start_index + proc_max]):
#             res = open(f'dock6/runs/{i}/flex.out').read().split('Secondary Score')[1]
#             smile_to_data[smile] = {'total_score': float(res.split('Grid_Score:')[1].split()[0].strip()),
#                                     'vdw_score': float(res.split('Grid_vdw_energy:')[1].split()[0].strip()),
#                                     'es_score': float(res.split('Grid_es_energy:')[1].split()[0].strip()),
#                                     'internal_energy': float(res.split('Internal_energy_repulsive:')[1].split()[0].strip())}
#             subprocess.call(f'obabel dock6/runs/{i}/dock1_secondary_scored.mol2 -O dock6/runs/{i}/out.pdb', shell=True)
#             atoms = []
#             for line in open(f'dock6/runs/{i}/out.pdb', 'r'):
#                 if line.startswith('ATOM'):
#                     _, _, type, _, _, x, y, z, vdw, _, _, _ = line.split()
#                     atoms.append((type, float(x), float(y), float(z)))
#             smile_to_data[smile]['atom_coords'] = atoms
#     return smile_to_data

# def dock6_amber(smiles):
#     smile_to_data = {}
#     for start_index in range(0, len(smiles), proc_max):
#         subprocess.run('rm -rf dock6/runs', shell=True, stderr=subprocess.DEVNULL)
#         os.mkdir('dock6/runs')
#         ps = []
#         for i, smile in enumerate(smiles[start_index:start_index + proc_max]):
#             dir = f'dock6/runs/{i}'
#             os.mkdir(dir)
#             ps.append(subprocess.Popen(f'obabel -:"{smile}" -O {dir}/ligand.mol2 -h --partialcharge mmff94 --gen3d', shell=True))
#         for p in ps:
#             p.wait()
#         open('run_dock6_amber_modified', 'w').write(open('run_dock6_amber', 'r').read().replace('_MAX_FILE_', str(len(smiles[start_index:start_index + proc_max]) - 1)))
#         subprocess.run(f'./run_dock6_amber_modified')
#         for i, smile in enumerate(smiles[start_index:start_index + proc_max]):
#             res = open(f'dock6/runs/{i}/amber.out').read().split('Conformations:')[1]
#             smile_to_data[smile] = {'amber_score': float(res.split('Amber_Score:')[1].split()[0].strip()),
#                                     'complex_energy': float(res.split('Amber_complex_energy:')[1].split()[0].strip()),
#                                     'receptor_energy': float(res.split('Amber_receptor_energy:')[1].split()[0].strip()),
#                                     'ligand_energy': float(res.split('Amber_ligand_energy:')[1].split()[0].strip())}
#             atoms = []
#             for line in open(f'dock6/runs/{i}/5uf0_noh.dock1_secondary_scored.1.final_pose.amber.pdb', 'r'):
#                 if line.startswith('ATOM'):
#                     _, _, type, _, _, x, y, z, vdw, _ = line.split()
#                     for char in '0123456789':
#                         type = type.replace(char, '')
#                     if type not in ['CL', 'BR']:
#                         type = type[0]
#                     atoms.append((type, float(x), float(y), float(z)))
#             smile_to_data[smile]['receptor_atom_coords'] = atoms
#             atoms = []
#             for line in open(f'dock6/runs/{i}/dock1_secondary_scored.1.final_pose.amber.pdb', 'r'):
#                 if line.startswith('ATOM'):
#                     _, _, type, _, _, x, y, z, vdw, _ = line.split()
#                     for char in '0123456789':
#                         type = type.replace(char, '')
#                     if type not in ['CL', 'BR']:
#                         type = type[0]
#                     atoms.append((type, float(x), float(y), float(z)))
#             smile_to_data[smile]['ligand_atom_coords'] = atoms
#     return smile_to_data

def surrogate(smiles):
    smile_to_data = {}
    x = smiles_to_morgan(smiles)
    out = surrogate_model.predict(x)
    for i, smile in enumerate(smiles):
        smile_to_data[smile] = {'prediction': out[i]}
    return smile_to_data