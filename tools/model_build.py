#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys, os
from pathlib import Path
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import time

current_file = Path(__file__).resolve()
project_root = current_file.parent.parent

if str(project_root) not in sys.path:
	sys.path.append(str(project_root))

from tools.io import STRSUMMARY
from parameters.load_params import readParam

HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))

def build_system(data, Mtype='DNA'):
	"""
		model construction
	"""

	# section of topology
	topology = app.Topology()
	info = STRSUMMARY(data, Mtype)
	
	atoms = []

	print('Starting build system: topology part')
	s_time = time.perf_counter()
	for chain in info['chains']:
		chain_tmp = topology.addChain() # add topology for current chain
		_dataframe_tmp = data[data['chain']==chain] # extract data belongs to current chain
		_Mtype = _dataframe_tmp['Mtype'].iloc[0] # since for this chain, it should belong to one Mtype
		_names = _dataframe_tmp['name'].to_dict()
		_indices = _dataframe_tmp.index.tolist()
		_bondlist = _dataframe_tmp['bondlist'].tolist()
		_atypes = _dataframe_tmp['ATYPE'].to_dict()

		if _Mtype == 'DNA':
			for _index in _indices:
				if _names[_index] == 'PO4':
					residue_tmp = topology.addResidue(name=_atypes[_index], chain = chain_tmp)
					atom_tmp    = topology.addAtom(name=_names[_index], element=None, residue=residue_tmp)
					atoms.append(atom_tmp)
				else:
					atom_tmp    = topology.addAtom(name=_names[_index], element=None, residue=residue_tmp)
					atoms.append(atom_tmp)
		elif _Mtype == 'protein':
			for _index in _indices:
				residue_tmp = topology.addResidue(name=_atypes[_index], chain=chain_tmp)
				atom_tmp    = topology.addAtom(name=_names[_index], element=None, residue=residue_tmp)
				atoms.append(atom_tmp)
		elif _Mtype == 'organic':
			pass
		else:
			raise('In building the system, however the molecule type {} is not implemented yet'.format(Mtype))

#		for i_x, row_x in data.iterrows():
#			Mtype = row_x['Mtype']
#			if row_x['chain'] == chain:
#				if Mtype == 'DNA':
#					#if i_x % 2 == 0: # probably should ancher the atomtype / name
#					if row_x['name'] == 'PO4':
#						residue_tmp = topology.addResidue(name=row_x['ATYPE'], chain=chain_tmp)
#						atom_tmp    = topology.addAtom(name=row_x['name'], element=None, residue=residue_tmp)
#						atoms.append(atom_tmp)
#					else:
#						atom_tmp    = topology.addAtom(name=row_x['name'], element=None, residue=residue_tmp)
#						atoms.append(atom_tmp)
#				elif Mtype == 'protein':
#					residue_tmp = topology.addResidue(name=row_x['ATYPE'], chain=chain_tmp)
#					atom_tmp    = topology.addAtom(name=row_x['name'], element=None, residue=residue_tmp)
#					atoms.append(atom_tmp)
#				elif Mtype == 'organic':
#					pass
#				else:
#					raise('In building the system, however the molecule type {} is not implemented yet'.format(Mtype))
	e_time = time.perf_counter()
	print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	print('Starting build system: bond part')
	s_time = time.perf_counter()
	for i, row in data.iterrows():
		if len(row['bondlist']):
			for nxt in row['bondlist']:
				topology.addBond(atoms[i], atoms[nxt])
	e_time = time.perf_counter()
	print('Elasping: {:.3f} seconds'.format(e_time - s_time))
	

	print('Starting build system: assign MASS and CHARGE')
	s_time = time.perf_counter()
	# section of system
	system = mm.System()
	for i, row in data.iterrows():
		if row['ATYPE'] in ['ADE', 'THY', 'CYT', 'GUA']:
			mass_tmp = HPS_DNA.loc[HPS_DNA['atomtype'] == row['name'],'mass'].values[0]
		else:
			mass_tmp = HPS_PROT.loc[HPS_PROT['atomtype']==row['ATYPE'],'mass'].values[0]
	
		#print(mass_tmp)
		system.addParticle(mass_tmp * unit.amu)

	e_time = time.perf_counter()
	print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	return topology, system
