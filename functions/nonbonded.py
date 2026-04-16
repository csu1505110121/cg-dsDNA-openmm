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

HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))


def stack_hb_exclude_DNA(dataframe, Mtype='DNA', use_pbc=True):
	
	# stack definition
	stack_expression = 'step(rstack -r) * (delta * (5* (sigma/r)^12 - 6*(sigma/r)^10))'
	stack_force      = mm.CustomBondForce(stack_expression)
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
	excluded_force      = mm.CustomBondForce(excluded_expression)
	excluded_force.addGlobalParameter('epsilon', 4*unit.kilocalorie_per_mole)
	excluded_force.addGlobalParameter('sigma', 0.55*unit.nanometer)

	# additional contribution from vdW
	vdw_expression = """
		step(r_vdw_cut - r)*step(r_vdw-r) * (4*epsilon*((sigma/r)^12 - (sigma/r)^6) + (1-lambda)*epsilon)
	"""
	vdw_force           = mm.CustomBondForce(vdw_expression)
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

	# ---- brute force style ---- #
#	for i, row_i in dataframe.iterrows():
#		Mtype_i = row_i['Mtype']
#		for j, row_j in dataframe.iterrows():
#			Mtype_j = row_j['Mtype']
#			if Mtype_i == 'DNA' and Mtype_j == 'DNA': # these force only apply to DNA
#				if i > j and not any(set([i,j]) == set(bond) for bond in bonds):
#					if Counter([row_i['name'],row_j['name']]) == Counter(['PO4','ADE']):
#						vdw_force.addBond(i, j)
#					elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','GUA']):
#						vdw_force.addBond(i, j)
#					elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','CYT']):
#						vdw_force.addBond(i, j)
#					elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','THY']):
#						vdw_force.addBond(i, j)
#					elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','PO4']):
#						vdw_force.addBond(i, j)
						
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
#
#	for i in dna_indices:
#		name_i = names[i]
#		chain_i = chains[i]
#
#		j2 = i-2
#		if j2 in names and chains[j2] == chain_i:
#			name_j2 = names[j2]
#			if name_i != 'PO4' and name_j2 != 'PO4':
#				pair = frozenset([name_i, name_j2])
#				if pair in stack_lookup:
#					eps_ = stack_lookup[pair]
#					stack_force.addBond(j2, i, [0.36*unit.nanometer, eps_*unit.kilocalorie_per_mole])
#		
#		j4 = i-4
#		if j4 in names and chains[j4] == chain_i:
#			if name_i != 'PO4' and names[j4] != 'PO4':
#				if name_i in target_names and names[j4] in target_names:
#					excluded_force.addBond(j4,i)
#
#	
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
				

###		name_i, name_j = names[i], names[j]
###		if name_i == 'PO4' or name_j == 'PO4':
###			continue
###
###		pair = frozenset([name_i, name_j])
###		if pair == frozenset(['ADE','THY']):
###			hb_force.addBond(j, i, [0.6*unit.nanometer, 2.7*unit.kilocalorie_per_mole])
###		elif pair == frozenset(['CYT', 'GUA']):
###			hb_force.addBond(j, i, [0.55*unit.nanometer, 3.3*unit.kilocalorie_per_mole])
###		# for exclude force
###		elif pair != frozenset(['ADE', 'THY']) and pair != frozenset(['CYT', 'GUA']):
###			if name_i in target_names and name_j in target_names:
###				excluded_force.addBond(j, i)
				

####  ---- brute force style ----
###	for i, row_i in dataframe.iterrows():
###		for j, row_j in dataframe.iterrows():
###			if i > j  and row_i['name'] != 'PO4' and row_j['name'] != 'PO4':
###				if row_i['chain'] == row_j['chain'] and abs(i-j) == 2: # make sure they are in the same chain and are adjacent to each base
###					if row_i['name'] == 'ADE' and row_j['name'] == 'ADE':
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 9.7*unit.kilocalorie_per_mole])
###					if row_i['name'] == 'THY' and row_j['name'] == 'THY':
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 7.7*unit.kilocalorie_per_mole])
###					if row_i['name'] == 'CYT' and row_j['name'] == 'CYT':
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 6.2*unit.kilocalorie_per_mole])
###					if row_i['name'] == 'GUA' and row_j['name'] == 'GUA':
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 11.0*unit.kilocalorie_per_mole])
###					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['ADE', 'THY']):
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 8.6*unit.kilocalorie_per_mole])
###					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['ADE', 'CYT']):
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 7.8*unit.kilocalorie_per_mole])
###					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['ADE', 'GUA']):
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 10.3*unit.kilocalorie_per_mole])
###					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['THY', 'CYT']):
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 6.9*unit.kilocalorie_per_mole])
###					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['THY', 'GUA']):
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 9.2*unit.kilocalorie_per_mole])
###					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['CYT', 'GUA']):
###						stack_force.addBond(i, j, [0.36*unit.nanometer, 8.3*unit.kilocalorie_per_mole])
###				elif row_i['chain'] == row_j['chain'] and abs(i-j) == 4:
###					if row_i['name'] in HPS_DNA['atomtype'].values and row_j['name'] in HPS_DNA['atomtype'].values: # just make sure for mixture: protein|DNA, protein will not contribute to exclude
###						excluded_force.addBond(i, j)
###				elif Counter([row_i['name'],row_j['name']]) != Counter(['ADE','THY']) and Counter([row_i['name'],row_j['name']]) != Counter(['CYT','GUA']):
###					if row_i['name'] in HPS_DNA['atomtype'].values and row_j['name'] in HPS_DNA['atomtype'].values: # just make sure for mixture: protein|DNA, protein will not contribute to exclude
###						excluded_force.addBond(i, j)
###				else:
###					if Counter([row_i['name'],row_j['name']]) == Counter(['ADE','THY']):
###						hb_force.addBond(i, j, [0.6*unit.nanometer,2.7*unit.kilocalorie_per_mole])
###					elif Counter([row_i['name'],row_j['name']]) == Counter(['CYT','GUA']):
###						hb_force.addBond(i, j, [0.55*unit.nanometer,3.3*unit.kilocalorie_per_mole])
						
	stack_force.setUsesPeriodicBoundaryConditions(use_pbc)
	stack_force.setForceGroup(3)

	hb_force.setUsesPeriodicBoundaryConditions(use_pbc)
	hb_force.setForceGroup(4)

	excluded_force.setUsesPeriodicBoundaryConditions(use_pbc)
	excluded_force.setForceGroup(5)
	
	return stack_force, hb_force, excluded_force, vdw_force


def vdw_12_6(dataframe, Mtype='protein', use_pbc=True, force_group=6):
	"""
		U_vdw = step(r_cut - r) * (4 * epsilon * [(\\frac{\\sigma}{r})^12 - (\\frac{\\sigma}{r})^6] + epsilon* (1-\\lambda)) + \
				step(r - r_cut) * (4 * epsilon * [(\\frac{\\sigma}{r})^12 - (\\frac{\\sigma}{r})^6])
	"""

	vdw_expression = """
		step(r_vdw_cut - r)* (lj + (1-lambda)*epsilon) + step(r-r_vdw_cut)*(lambda*lj); 
											lj = 4*epsilon*((sigma/r)^12 - (sigma/r)^6); 
											r_vdw_cut = 2^(1/6)*sigma; 
											epsilon = (epsilon1 * epsilon2)^0.5;
											lambda = 0.5*(lambda1 + lambda2); 
											sigma=0.5*(sigma1+sigma2)
	"""

	vdw_force           = mm.CustomNonbondedForce(vdw_expression)
	vdw_force.addPerParticleParameter('epsilon')
	vdw_force.addPerParticleParameter('sigma')
	vdw_force.addPerParticleParameter('lambda')


	sigma_dna_map = HPS_DNA.set_index('atomtype')['sigma'].to_dict()
	lambd_dna_map = HPS_DNA.set_index('atomtype')['lambda'].to_dict()

	sigma_prot_map = HPS_PROT.set_index('atomtype')['sigma'].to_dict()
	lambd_prot_map = HPS_PROT.set_index('atomtype')['lambda'].to_dict()

	# store all bonds info including both 'DNA' and 'protein'
	group_protein = dataframe[dataframe['Mtype']=='protein'].index.tolist() # store index for protein group
	group_DNA     = dataframe[dataframe['Mtype']=='DNA'].index.tolist() # store index for DNA group
	_bondlists    = dataframe['bondlist'].tolist()
	_indices      = dataframe.index.tolist()
	_Mtypes       = dataframe['Mtype'].tolist()
	_names        = dataframe['name'].tolist()
	_atypes       = dataframe['ATYPE'].tolist()

	bonds = []
	for _indice, _bondlist in zip(_indices, _bondlists):
		if len(_bondlist):
			for nxt in _bondlist:
				bonds.append((_indice, nxt))
			

####  old script

#	for x, row_x in dataframe.iterrows():
#		if row_x['Mtype'] == 'protein':
#			group_protein.append(x)
#		elif row_x['Mtype'] == 'DNA':
#			group_DNA.append(x)
#		if len(row_x['bondlist']):
#			for nxt in row_x['bondlist']:
#				bonds.append((x, nxt))

	max_radii = 0
	for _mtype, _name, _atype in zip(_Mtypes, _names, _atypes):
		if _mtype == 'DNA':
			eps = 0.2
			radii = sigma_dna_map.get(_name)
			lba   = lambd_dna_map.get(_name)
		elif _mtype == 'protein':
			eps = 0.2
			radii = sigma_prot_map.get(_atype)
			lba   = lambd_prot_map.get(_atype)
		if radii > max_radii:
			max_radii = radii

		vdw_force.addParticle([eps* unit.kilocalorie_per_mole, radii* unit.nanometer, lba])



	## old script
###	# assign parameters
###	max_radii = 0
###	for i, row in dataframe.iterrows():
###		Mtype = row['Mtype']
###		if Mtype == 'DNA':
###			eps   = 0.2
###			radii = HPS_DNA.loc[HPS_DNA['atomtype']==row['name'],'sigma'].values[0]
###			lba = HPS_DNA.loc[HPS_DNA['atomtype']==row['name'],'lambda'].values[0]
###		elif Mtype == 'protein':
###			eps  = 0.2
###			radii = HPS_PROT.loc[HPS_PROT['atomtype']==row['ATYPE'],'sigma'].values[0]
###			lba   = HPS_PROT.loc[HPS_PROT['atomtype']==row['ATYPE'],'lambda'].values[0]
###		if radii > max_radii:
###			max_radii = radii
###		# three parameters are required:
###		# epsilon
###		# sigma
###		# lambda
###		vdw_force.addParticle([eps* unit.kilocalorie_per_mole, radii* unit.nanometer, lba])
		
	if len(group_protein) !=0 and len(group_DNA) !=0:
		vdw_force.addInteractionGroup(group_protein, group_protein)
		vdw_force.addInteractionGroup(group_protein, group_DNA)
	elif len(group_protein) and len(group_DNA) ==0:
		vdw_force.addInteractionGroup(group_protein, group_protein)


	vdw_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
	vdw_force.createExclusionsFromBonds(bonds,1)
	vdw_force.setCutoffDistance(4*max_radii * unit.nanometer)
	vdw_force.setForceGroup(force_group)

	return vdw_force


def elec(dataframe, Mtype='DNA', dielectric=80, use_pbc=True, force_group=7):
	"""
		U_electrostatics = \\frac{q_i * q_j}{4*\\pi*\\epsilon_0 * \\epsilon_r * r_{ij}} \\exp^{-r_ij/\\kappa}
	"""

	NA=unit.AVOGADRO_CONSTANT_NA
	EC=1.602176634e-19*unit.coulomb
	VEP=8.8541878128e-12*unit.farad/unit.meter
	alpha = NA*EC**2 / (4*np.pi*VEP)
	ldby=1*unit.nanometer
	alpha_value=alpha.value_in_unit(unit.kilojoule_per_mole*unit.nanometer)
	ldby_value = ldby.value_in_unit(unit.nanometer)
	elec_expression = 'q1*q2*alpha_value*(exp(-r/ldby_value))/r/dielectric_water'
	elec_force = mm.CustomNonbondedForce(elec_expression)

	elec_force.addGlobalParameter('alpha_value',alpha)
	elec_force.addGlobalParameter('ldby_value',ldby)
	elec_force.addGlobalParameter('dielectric_water',dielectric)

	elec_force.addPerParticleParameter('q')


	chg_dna_map = HPS_DNA.set_index('atomtype')['charge'].to_dict()
	chg_prot_map = HPS_PROT.set_index('atomtype')['charge'].to_dict()



	_indices      = dataframe.index.tolist()
	_Mtypes       = dataframe['Mtype'].tolist()
	_names        = dataframe['name'].tolist()
	_atypes       = dataframe['ATYPE'].tolist()
	_bondlists    = dataframe['bondlist'].tolist()

	for _mtype, _name, _atype in zip(_Mtypes, _names, _atypes):
		if _mtype == 'DNA':
			chg = chg_dna_map.get(_name)
			elec_force.addParticle([chg])
		elif _mtype == 'protein':
			chg = chg_prot_map.get(_atype)
			elec_force.addParticle([chg])
		else:
			raise('In applying elec force, however the molecule type {} is not implemented yet'.format(Mtype))
			
#### old script
#	for i, row in dataframe.iterrows():
#		Mtype = row['Mtype']
#		if Mtype == 'DNA':
#			chg = HPS_DNA.loc[HPS_DNA['atomtype']==row['name'],'charge'].values[0]
#			elec_force.addParticle([chg])
#		elif Mtype == 'protein':
#			chg = HPS_PROT.loc[HPS_PROT['atomtype']==row['ATYPE'],'charge'].values[0]
#			elec_force.addParticle([chg])
#		else:
#			raise('In applying elec force, however the molecule type {} is not implemented yet'.format(Mtype))
			
	bonds = []
	for _indice, _bondlist in zip(_indices, _bondlists):
		if len(_bondlist):
			for nxt in _bondlist:
				bonds.append((_indice, nxt))

#	# store all bonds info
#	bonds = []
#	for x, row_x in dataframe.iterrows():
#		if len(row_x['bondlist']):
#			for nxt in row_x['bondlist']:
#				bonds.append((x, nxt))

	# should pay attention to this line
	#elec_force.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()],1)
	elec_force.createExclusionsFromBonds(bonds,1)

	elec_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
	elec_force.setCutoffDistance(3.5*unit.nanometer)
	elec_force.setForceGroup(force_group)

	return elec_force

if __name__ == '__main__':
	#print(current_file.parent.parent)
	print(HPS_PROT)

