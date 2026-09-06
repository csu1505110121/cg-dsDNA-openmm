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



def build_system(data:pd.DataFrame, info_conf:dict):
	"""
		model construction
	"""
	HPS_PROT = readParam(info_conf['fn_hps_prot'])
	HPS_DNA = readParam(info_conf['fn_hps_DNA'])
	HPS_NUCLEOTIDE = readParam(info_conf['fn_hps_nucleotide'])	

	# section of topology
	topology = app.Topology()
	info = STRSUMMARY(data)
	
	atoms = []

	print('Starting build system: topology part')
	s_time = time.perf_counter()
	for chain in info['chains']:
		chain_tmp = topology.addChain() # add topology for current chain
		_dataframe_tmp = data[data['chain']==chain] # extract data belongs to current chain
		_Mtype = _dataframe_tmp['Mtype'].iloc[0] # since for this chain, it should belong to one Mtype
		_names = _dataframe_tmp['name'].to_dict()
		_indices = _dataframe_tmp.index.tolist()
		#_bondlist = _dataframe_tmp['bondlist'].tolist()
		_atypes = _dataframe_tmp['ATYPE'].to_dict()

		#print(_names, _atypes)

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
		elif _Mtype == 'nucleotide':
			for _index in _indices:
				residue_tmp = topology.addResidue(name=_atypes[_index], chain=chain_tmp)
				atom_tmp    = topology.addAtom(name=_names[_index], element=None, residue=residue_tmp)
				atoms.append(atom_tmp)
		else:
			raise('In building the system, however the molecule type {} is not implemented yet'.format(_Mtype))

	e_time = time.perf_counter()
	print('Elasped: {:.3f} seconds'.format(e_time - s_time))

	print('Starting build system: bond part')
	s_time = time.perf_counter()
	_bondlists = data['bondlist'].tolist()
	_indices  = data.index.tolist()
	_atypes   = data['ATYPE'].tolist()
	_names    = data['name'].tolist()

	#print(_atypes)

	for _index, _bondlist in zip(_indices, _bondlists):
		if len(_bondlist):
			for nxt in _bondlist:
				topology.addBond(atoms[_index], atoms[nxt])
	e_time = time.perf_counter()
	print('Elasped: {:.3f} seconds'.format(e_time - s_time))


	mass_dna_map = HPS_DNA.set_index('atomtype')['mass'].to_dict()
	mass_prot_map = HPS_PROT.set_index('atomtype')['mass'].to_dict()
	mass_nucleotide_map = HPS_NUCLEOTIDE.set_index('atomtype')['mass'].to_dict()

	#print(mass_nucleotide_map)

	print('Starting build system: assign MASS part')
	s_time = time.perf_counter()
	# section of system
	system = mm.System()
	for _index, _atype, _name in zip(_indices, _atypes, _names):
		if _atype in ['ADE', 'THY', 'CYT', 'GUA']:
			mass_tmp = mass_dna_map.get(_atype)
		elif _atype in ['ATP', 'ADP', 'AMP', 'CTP', 'CDP', 'CMP', 'GTP', 'GDP', 'GMP', 'TTP', 'TDP', 'TMP']: # need another identifier for NUC
			mass_tmp = mass_nucleotide_map.get(_name)
		else:
			mass_tmp = mass_prot_map.get(_atype)
	
		#print(mass_tmp)
		system.addParticle(mass_tmp * unit.amu)
	e_time = time.perf_counter()
	print('Elasped: {:.3f} seconds'.format(e_time - s_time))


	return topology, system
