#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from collections import Counter
import sys,os
import itertools
from pathlib import Path
import openmm as mm
import openmm.app as app
import openmm.unit as unit

current_file = Path(__file__).resolve()
project_root = current_file.parent.parent

if str(project_root) not in sys.path:
	sys.path.append(str(project_root))

from parameters.load_params import readParam

#HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
#HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))
#HPS_NUCLEOTIDE = readParam(os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv'))


def stack_hb_exclude_DNA(dataframe, use_pbc=True):
	
	# stack definition
	stack_expression = 'step(rstack -r) * (delta * (5* (sigma/r)^12 - 6*(sigma/r)^10))'
	stack_force	  = mm.CustomBondForce(stack_expression)
	stack_force.addGlobalParameter('rstack', 0.62*unit.nanometer)

	stack_force.addPerBondParameter('sigma')
	stack_force.addPerBondParameter('delta')

	# hb definition
	hb_expression = """
				step(r_hb - r)*(delta_hb*(5*(sigma_hb/r)^12 - 6*(sigma_hb/r)^10))
	"""
	hb_force = mm.CustomBondForce(hb_expression)
	hb_force.addGlobalParameter('r_hb', 0.95*unit.nanometer)

	hb_force.addPerBondParameter('sigma_hb')
	hb_force.addPerBondParameter('delta_hb')
	
	# excluded definition
	excluded_expression = """
			step(rcut-r)*(4*epsilon*((sigma/r)^12 - (sigma/r)^6) + epsilon); rcut=2^(1/6)*sigma
	"""
	excluded_force	  = mm.CustomBondForce(excluded_expression)
	excluded_force.addGlobalParameter('epsilon', 4*unit.kilocalorie_per_mole)
	excluded_force.addGlobalParameter('sigma', 0.55*unit.nanometer)

	# additional contribution from vdW
	vdw_expression = """
		step(r_vdw_cut - r)*step(r_vdw-r) * (4*epsilon*((sigma/r)^12 - (sigma/r)^6) + (1-lambda)*epsilon)
	"""
	vdw_force		   = mm.CustomBondForce(vdw_expression)
	vdw_force.addGlobalParameter('lambda',0)
	vdw_force.addGlobalParameter('r_vdw',0.55*2**(1/6)*unit.nanometer)
	vdw_force.addGlobalParameter('r_vdw_cut',0.6174*unit.nanometer)

	vdw_force.addGlobalParameter('epsilon',4*unit.kilocalorie_per_mole)
	vdw_force.addGlobalParameter('sigma',0.55*unit.nanometer)

	
#	# store all bonds info
#	bonds = []
#	for x, row_x in dataframe.iterrows():
#		Mtype = row_x['Mtype']
#		if Mtype == 'DNA':
#			if len(row_x['bondlist']):
#				for nxt in row_x['bondlist']:
#					bonds.append([x, nxt])

	# store Mtype == DNA into another dataframe
	dna_df = dataframe[dataframe['Mtype'] == 'DNA']
	names = dna_df['name'].to_dict()
	chains = dna_df['chain'].to_dict()
	dna_indices = dna_df.index.tolist()
	dna_bondlist = dna_df['bondlist'].tolist()

	# store dsDNA bonds info
	bonds = []
	for x, bond_x in zip(dna_indices, dna_bondlist):
		if len(bond_x):
			for nxt in bond_x:
				bonds.append([x, nxt])


	existing_bonds = {frozenset(b) for b in bonds}
	target_names = {'ADE','GUA','CYT','THY','PO4'}

	for i,j in itertools.combinations(dna_indices, 2):
		if frozenset([i,j]) in existing_bonds:
			continue

		name_i, name_j = names[i], names[j]

		if (name_i == 'PO4' and name_j in target_names) or (name_j == 'PO4' and name_i in target_names):
			vdw_force.addBond(i,j)

						
	vdw_force.setUsesPeriodicBoundaryConditions(use_pbc)
	vdw_force.setForceGroup(6)


	stack_lookup = {
		frozenset(['ADE', 'ADE']): 9.7,
		frozenset(['THY', 'THY']): 7.7,
		frozenset(['CYT', 'CYT']): 6.2,
		frozenset(['GUA', 'GUA']): 11.0,
		frozenset(['ADE', 'THY']): 8.6,
		frozenset(['ADE', 'CYT']): 7.8,
		frozenset(['ADE', 'GUA']): 10.3,
		frozenset(['THY', 'CYT']): 6.9,
		frozenset(['THY', 'GUA']): 9.2,
		frozenset(['CYT', 'GUA']): 8.3
	}

	
	for i,j in itertools.combinations(dna_indices, 2):
		diff = abs(i-j)
		name_i, name_j = names[i], names[j]
		chain_i, chain_j = chains[i], chains[j]
		if name_i != 'PO4' and name_j != 'PO4':
			pair = frozenset([name_i, name_j])
			if diff == 2 and chain_i == chain_j:
				if pair in stack_lookup:
					eps_ = stack_lookup[pair]
					stack_force.addBond(i,j, [0.36*unit.nanometer, eps_*unit.kilocalorie_per_mole])
		
			elif diff == 4 and chain_i == chain_j:
				excluded_force.addBond(i,j)
	
			elif pair != frozenset(['ADE','THY']) and pair != frozenset(['CYT','GUA']):
				excluded_force.addBond(i,j)
			else:
				if pair == frozenset(['ADE','THY']):
					hb_force.addBond(i, j, [0.6*unit.nanometer, 2.7*unit.kilocalorie_per_mole])
				elif pair == frozenset(['CYT', 'GUA']):
					hb_force.addBond(i, j, [0.55*unit.nanometer, 3.3*unit.kilocalorie_per_mole])
				

						
	stack_force.setUsesPeriodicBoundaryConditions(use_pbc)
	stack_force.setForceGroup(3)

	hb_force.setUsesPeriodicBoundaryConditions(use_pbc)
	hb_force.setForceGroup(4)

	excluded_force.setUsesPeriodicBoundaryConditions(use_pbc)
	excluded_force.setForceGroup(5)
	
	return stack_force, hb_force, excluded_force, vdw_force


def vdw_12_6(dataframe, info_conf, use_pbc=True, force_group=6):
	"""
		U_vdw = step(r_cut - r) * (4 * epsilon * [(\\frac{\\sigma}{r})^12 - (\\frac{\\sigma}{r})^6] + epsilon* (1-\\lambda)) + \
				step(r - r_cut) * (4 * epsilon * [(\\frac{\\sigma}{r})^12 - (\\frac{\\sigma}{r})^6])
	"""

	if 'fn_hps_prot' in info_conf:
		file = Path(info_conf['fn_hps_prot'])
		if file.exists():
			HPS_PROT = readParam(info_conf['fn_hps_prot'])
		else:
			HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
	else:
		HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
	if 'fn_hps_DNA' in info_conf:
		file = Path(info_conf['fn_hps_DNA'])
		if file.exists():
			HPS_DNA = readParam(info_conf['fn_hps_DNA'])
		else:
			HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))
	else:
		HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))
	if 'fn_hps_nucleotide' in info_conf:
		file = Path(info_conf['fn_hps_nucleotide'])
		if file.exists():
			HPS_NUCLEOTIDE = readParam(info_conf['fn_hps_nucleotide'])
		else:
			HPS_NUCLEOTIDE = readParam(os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv'))
	else:
		HPS_NUCLEOTIDE = readParam(os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv'))



	vdw_expression = """
		mask_val * step(r_vdw_cut - r)* (lj + (1-lambda)*epsilon) + mask_val * step(r-r_vdw_cut)*(lambda*lj); 
											lj = 4*epsilon*((sigma/r)^12 - (sigma/r)^6); 
											r_vdw_cut = 2^(1/6)*sigma; 
											epsilon = (epsilon1 * epsilon2)^0.5;
											lambda = 0.5*(lambda1 + lambda2); 
											sigma=0.5*(sigma1+sigma2);
											mask_val = step(1.5 - (is_dna1 + is_dna2))
	"""

	vdw_force		   = mm.CustomNonbondedForce(vdw_expression)
	vdw_force.addPerParticleParameter('epsilon')
	vdw_force.addPerParticleParameter('sigma')
	vdw_force.addPerParticleParameter('lambda')
	vdw_force.addPerParticleParameter('is_dna')

	sigma_nucleotide_map = HPS_NUCLEOTIDE.set_index('atomtype')['sigma'].to_dict()
	lambd_nucleotide_map = HPS_NUCLEOTIDE.set_index('atomtype')['lambda'].to_dict()
	epsilon_nucleotide_map = HPS_NUCLEOTIDE.set_index('atomtype')['epsilon'].to_dict()

	sigma_dna_map = HPS_DNA.set_index('atomtype')['sigma'].to_dict()
	lambd_dna_map = HPS_DNA.set_index('atomtype')['lambda'].to_dict()
	epsilon_dna_map = HPS_DNA.set_index('atomtype')['epsilon'].to_dict()

	sigma_prot_map = HPS_PROT.set_index('atomtype')['sigma'].to_dict()
	lambd_prot_map = HPS_PROT.set_index('atomtype')['lambda'].to_dict()
	epsilon_prot_map = HPS_PROT.set_index('atomtype')['epsilon'].to_dict()

	# store all bonds info including both 'DNA' and 'protein'
	group_protein = dataframe[dataframe['Mtype']=='protein'].index.tolist() # store index for protein group
	group_DNA	 = dataframe[dataframe['Mtype']=='DNA'].index.tolist() # store index for DNA group
	_bondlists	= dataframe['bondlist'].tolist()
	_indices	  = dataframe.index.tolist()
	_Mtypes	   = dataframe['Mtype'].tolist()
	_names		= dataframe['name'].tolist()
	_atypes	   = dataframe['ATYPE'].tolist()

	bonds = []
	for _indice, _bondlist in zip(_indices, _bondlists):
		if len(_bondlist):
			for nxt in _bondlist:
				bonds.append((_indice, nxt))
			


	max_radii = 0
	for _mtype, _name, _atype in zip(_Mtypes, _names, _atypes):
		if _mtype == 'DNA':
			eps = epsilon_dna_map.get(_name)
			radii = sigma_dna_map.get(_name)
			lba   = lambd_dna_map.get(_name)
			is_dna = 1
		elif _mtype == 'protein':
			eps = epsilon_prot_map.get(_atype)
			radii = sigma_prot_map.get(_atype)
			lba   = lambd_prot_map.get(_atype)
			is_dna = 0
		elif _mtype == 'nucleotide':
			eps = epsilon_nucleotide_map.get(_name)
			radii = sigma_nucleotide_map.get(_name)
			lba   = lambd_nucleotide_map.get(_name)
			is_dna = 0
		if radii > max_radii:
			max_radii = radii

		vdw_force.addParticle([eps* unit.kilocalorie_per_mole, radii* unit.nanometer, lba, is_dna])

	if use_pbc:
		vdw_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
	vdw_force.createExclusionsFromBonds(bonds,1)
	vdw_force.setCutoffDistance(4*max_radii * unit.nanometer)
	vdw_force.setForceGroup(force_group)

	print('CUTOFF distance for vdw force is set to {} nm'.format(4*max_radii))

	return vdw_force


def elec(dataframe, info_conf, dielectric=80, C_salt=1, use_pbc=True, force_group=7):
	"""
		U_electrostatics = \\frac{q_i * q_j}{4*\\pi*\\epsilon_0 * \\epsilon_r * r_{ij}} \\exp^{-r_ij/\\kappa}
	"""

	NA=unit.AVOGADRO_CONSTANT_NA
	EC=1.602176634e-19*unit.coulomb
	VEP=8.8541878128e-12*unit.farad/unit.meter
	alpha = NA*EC**2 / (4*np.pi*VEP)
	ldby=C_salt*unit.nanometer
	alpha_value=alpha.value_in_unit(unit.kilojoule_per_mole*unit.nanometer)
	ldby_value = ldby.value_in_unit(unit.nanometer)
	elec_expression = 'q1*q2*alpha_value*(exp(-r/ldby_value))/r/dielectric_water'
	elec_force = mm.CustomNonbondedForce(elec_expression)

	elec_force.addGlobalParameter('alpha_value',alpha_value)
	elec_force.addGlobalParameter('ldby_value',ldby_value)
	elec_force.addGlobalParameter('dielectric_water',dielectric)

	elec_force.addPerParticleParameter('q')


	if 'fn_hps_prot' in info_conf:
		file = Path(info_conf['fn_hps_prot'])
		if file.exists():
			HPS_PROT = readParam(info_conf['fn_hps_prot'])
		else:
			HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
	else:
		HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
	if 'fn_hps_DNA' in info_conf:
		file = Path(info_conf['fn_hps_DNA'])
		if file.exists():
			HPS_DNA = readParam(info_conf['fn_hps_DNA'])
		else:
			HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))
	else:
		HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))
	if 'fn_hps_nucleotide' in info_conf:
		file = Path(info_conf['fn_hps_nucleotide'])
		if file.exists():
			HPS_NUCLEOTIDE = readParam(info_conf['fn_hps_nucleotide'])
		else:
			HPS_NUCLEOTIDE = readParam(os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv'))
	else:
		HPS_NUCLEOTIDE = readParam(os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv'))


	chg_dna_map = HPS_DNA.set_index('atomtype')['charge'].to_dict()
	chg_prot_map = HPS_PROT.set_index('atomtype')['charge'].to_dict()
	chg_nucleotide_map = HPS_NUCLEOTIDE.set_index('atomtype')['charge'].to_dict()



	_indices	  = dataframe.index.tolist()
	_Mtypes	   = dataframe['Mtype'].tolist()
	_names		= dataframe['name'].tolist()
	_atypes	   = dataframe['ATYPE'].tolist()
	_bondlists	= dataframe['bondlist'].tolist()

	for _mtype, _name, _atype in zip(_Mtypes, _names, _atypes):
		if _mtype == 'DNA':
			chg = chg_dna_map.get(_name)
			elec_force.addParticle([chg])
		elif _mtype == 'protein':
			chg = chg_prot_map.get(_atype)
			elec_force.addParticle([chg])
		elif _mtype == 'nucleotide':
			chg = chg_nucleotide_map.get(_name)
			elec_force.addParticle([chg])
		else:
			raise('In applying elec force, however the molecule type {} is not implemented yet'.format(Mtype))
			
			
	bonds = []
	for _indice, _bondlist in zip(_indices, _bondlists):
		if len(_bondlist):
			for nxt in _bondlist:
				bonds.append((_indice, nxt))


	# should pay attention to this line
	#elec_force.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()],1)
	elec_force.createExclusionsFromBonds(bonds,1)

	if use_pbc:
		elec_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
	elec_force.setCutoffDistance(3.5*unit.nanometer)
	elec_force.setForceGroup(force_group)

	print('CUTOFF distance for elec force is set to {} nm'.format(3.5))	
	return elec_force

if __name__ == '__main__':
	#print(current_file.parent.parent)
	print(HPS_PROT)

