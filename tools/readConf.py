#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys, os
from pathlib import Path

current_file = Path(__file__).resolve()
project_root = current_file.parent.parent

if str(project_root) not in sys.path:
	sys.path.append(str(project_root))




def parseConf(fn,summary=False, ofn='mdp'):
	info = {
			# structure loading
			'xyz':[],
			'protein':[],
			'DNA':[],
			'nucleotide':[],
			'fn_protein':[],
			'fn_DNA':[],
			'fn_nucleotide':[],
			# md parameters
			'interval':[],
			'timestep':[],
			'temperature':[],
			'dielectric':[80],
			'cellVector1':[],
			'cellVector2':[],
			'cellVector3':[],
			'cellOrigin':[],
			'anisotropy':[False,False,False],
			# machine/platform specify
			'platform':['CUDA'],
			'prop':{'Threads':'1'},
			# output log file
			'log':'log.openmm',
			'trajfile':'trajtest.dcd',
			'rstfile': 'checkpoint.chk',
			# advanced parameters
			'enhanced_sampling':0,
			'plumed_conf':'',
			}
	with open(fn,'r') as f:
		while True:
			line = f.readline()
			if not line:
				break
			else:
				if '#' not in line and len(line.strip()):
					_record = line.split(';')[0]
					if len(_record.split()):
						record  = _record.split('=')
						args, vals = record[0].strip(), record[1].strip()
						try:
							if '.' in vals:
								vals = float(vals)
							else:
								vals = int(vals)
						except ValueError:
							pass
						info[args] = vals
	# convert to correct platform and corresponding properties
	if info['anisotropy']:
		if type(info['anisotropy']) == str:
			info['anisotropy'] = [x.strip().lower() == 'true' for x in info['anisotropy'].split(',')]

	if info['platform']	 in ['CUDA','HIP']:
		if 'DeviceIndex' in info.keys() and len(str(info['DeviceIndex'])):
			deviceidx = str(info['DeviceIndex'])
			info.pop('DeviceIndex')
		else:
			deviceidx = '0' # a default value
		if 'Precision' in info.keys() and len(info['Precision']):
			precision = info['Precision']
			info.pop('Precision')
		else:
			precision = 'mixed'
		info['prop'] = {'DeviceIndex':deviceidx, 'Precision': precision}

	if info['platform'] == 'CPU':
		if 'Threads' in info.keys() and len(str(info['Threads'])):
			threads = str(info['Threads'])
			info.pop('Threads')
		else:
			threads = '1' # set default value to 1

		info['prop'] = {'Threads':threads}

	if len(info['cellVector1']) and len(info['cellVector2']) and len(info['cellVector3']):
		v1 = [float(x) for x in info['cellVector1'].split(',')]
		v2 = [float(x) for x in info['cellVector2'].split(',')]
		v3 = [float(x) for x in info['cellVector3'].split(',')]
		info['pbc'] = ((v1[0],v1[1],v1[2]),(v2[0],v2[1],v2[2]),(v3[0],v3[1],v3[2]))
		info.pop('cellVector1')
		info.pop('cellVector2')
		info.pop('cellVector3')

	## in case multiple protein chains are provided, need to convert them into a list of integers
	if info['protein']:
		if type(info['protein']) == str:
			info['protein'] = [int(x) for x in info['protein'].split(',')]
			info['fn_protein']	= [x.strip() for x in info['fn_protein'].split(',')]
		else:
			info['protein'] = [info['protein']]
			info['fn_protein'] = [info['fn_protein']]
	else:
		info['protein'] = [int(0)]

	if info['DNA']:
		if type(info['DNA']) == str:
			info['DNA'] = [int(x) for x in info['DNA'].split(',')]
			info['fn_DNA']	= [x.strip() for x in info['fn_DNA'].split(',')]
		else:
			info['DNA'] = [info['DNA']]
			info['fn_DNA'] = [info['fn_DNA']]
	else:
		info['DNA'] = [int(0)]

	if info['nucleotide']:
		if type(info['nucleotide']) == str:
			info['nucleotide'] = [int(x) for x in info['nucleotide'].split(',')]
			info['fn_nucleotide']	= [x.strip() for x in info['fn_nucleotide'].split(',')]
		else:
			info['nucleotide'] = [info['nucleotide']]
			info['fn_nucleotide'] = [info['fn_nucleotide']]
	else:
		info['nucleotide'] = [int(0)]

	# read the default parameters or user-defined parameters in configure file

	if info['protein']:
		if 'fn_hps_prot' in info.keys() and len(info['fn_hps_prot']):
			file = Path(info['fn_hps_prot'])
			if not file.exists():
				print(f"File {info['fn_hps_prot']} does not exist and use the default parameters.")
				info['fn_hps_prot'] = os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv')
		else:
			info['fn_hps_prot'] = os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv')
	else:
		info['fn_hps_prot'] = os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv')


	if info['DNA']:
		if 'fn_hps_DNA' in info.keys() and len(info['fn_hps_DNA']):
			file = Path(info['fn_hps_DNA'])
			if not file.exists():
				print(f"File {info['fn_hps_DNA']} does not exist and use the default parameters.")
				info['fn_hps_DNA'] = os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv')
		else:
			info['fn_hps_DNA'] = os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv')
	else:
		info['fn_hps_DNA'] = os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv')


	if info['nucleotide']:
		if 'fn_hps_nucleotide' in info.keys() and len(info['fn_hps_nucleotide']):
			file = Path(info['fn_hps_nucleotide'])
			if not file.exists():
				print(f"File {info['fn_hps_nucleotide']} does not exist and use the default parameters.")
				info['fn_hps_nucleotide'] = os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv')
		else:
			info['fn_hps_nucleotide'] = os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv')

		if 'fn_nucleotide_stacking' in info.keys() and len(info['fn_nucleotide_stacking']):
			file = Path(info['fn_nucleotide_stacking'])
			if not file.exists():
				print(f"File {info['fn_nucleotide_stacking']} does not exist and use the default parameters.")
				info['fn_nucleotide_stacking'] = os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides_stacking.csv')
		else:
			info['fn_nucleotide_stacking'] = os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides_stacking.csv')
	else:
		info['fn_hps_nucleotide'] = os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv')
		info['fn_nucleotide_stacking'] = os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides_stacking.csv')



	if summary:
		with open(ofn, 'w') as f:
			for k,v in info.items():
				f.write('{:10s} \t = {}\n'.format(k,v))
		

	return info


if __name__ == '__main__':
	fn_template = os.path.join(project_root,'template/md.conf')
	info = parseConf(fn_template)
	print(info)
