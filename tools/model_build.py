#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys, os
from pathlib import Path
import openmm as mm
import openmm.app as app
import openmm.unit as unit

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

	for chain in info['chains']:
		chain_tmp = topology.addChain()
		for i_x, row_x in data.iterrows():
		#	if len(row_x['bondlist']):
		#		for nxt in row_x['bondlist']:
		#			topology.addBond(i_x, nxt)
			Mtype = row_x['Mtype']
			if row_x['chain'] == chain:
				if Mtype == 'DNA':
					#if i_x % 2 == 0: # probably should ancher the atomtype / name
					if row_x['name'] == 'PO4':
						residue_tmp = topology.addResidue(name=row_x['ATYPE'], chain=chain_tmp)
						atom_tmp    = topology.addAtom(name=row_x['name'], element=None, residue=residue_tmp)
						atoms.append(atom_tmp)
					else:
						atom_tmp    = topology.addAtom(name=row_x['name'], element=None, residue=residue_tmp)
						atoms.append(atom_tmp)
				elif Mtype == 'protein':
					residue_tmp = topology.addResidue(name=row_x['ATYPE'], chain=chain_tmp)
					atom_tmp    = topology.addAtom(name=row_x['name'], element=None, residue=residue_tmp)
					atoms.append(atom_tmp)
				elif Mtype == 'organic':
					pass
				else:
					raise('In building the system, however the molecule type {} is not implemented yet'.format(Mtype))

	for i, row in data.iterrows():
		if len(row['bondlist']):
			for nxt in row['bondlist']:
				topology.addBond(atoms[i], atoms[nxt])

	# section of system
	system = mm.System()
	for i, row in data.iterrows():
		if row['ATYPE'] in ['ADE', 'THY', 'CYT', 'GUA']:
			mass_tmp = HPS_DNA.loc[HPS_DNA['atomtype'] == row['name'],'mass'].values[0]
		else:
			mass_tmp = HPS_PROT.loc[HPS_PROT['atomtype']==row['ATYPE'],'mass'].values[0]
	
		#print(mass_tmp)
		system.addParticle(mass_tmp * unit.amu)

	return topology, system
