#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from collections import Counter
import sys,os
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

	
	# store all bonds info
	bonds = []
	for x, row_x in dataframe.iterrows():
		if len(row_x['bondlist']):
			for nxt in row_x['bondlist']:
				bonds.append([x, nxt])

	for i, row_i in dataframe.iterrows():
		for j, row_j in dataframe.iterrows():
			if i > j and not any(set([i,j]) == set(bond) for bond in bonds):
				if Counter([row_i['name'],row_j['name']]) == Counter(['PO4','ADE']):
					vdw_force.addBond(i, j)
				elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','GUA']):
					vdw_force.addBond(i, j)
				elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','CYT']):
					vdw_force.addBond(i, j)
				elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','THY']):
					vdw_force.addBond(i, j)
				elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','PO4']):
					vdw_force.addBond(i, j)
						
	vdw_force.setUsesPeriodicBoundaryConditions(use_pbc)
	vdw_force.setForceGroup(6)

	
	for i, row_i in dataframe.iterrows():
		for j, row_j in dataframe.iterrows():
			if i > j  and row_i['name'] != 'PO4' and row_j['name'] != 'PO4':
				if row_i['chain'] == row_j['chain'] and abs(i-j) == 2: # make sure they are in the same chain and are adjacent to each base
					if row_i['name'] == 'ADE' and row_j['name'] == 'ADE':
						stack_force.addBond(i, j, [0.36*unit.nanometer, 9.7*unit.kilocalorie_per_mole])
					if row_i['name'] == 'THY' and row_j['name'] == 'THY':
						stack_force.addBond(i, j, [0.36*unit.nanometer, 7.7*unit.kilocalorie_per_mole])
					if row_i['name'] == 'CYT' and row_j['name'] == 'CYT':
						stack_force.addBond(i, j, [0.36*unit.nanometer, 6.2*unit.kilocalorie_per_mole])
					if row_i['name'] == 'GUA' and row_j['name'] == 'GUA':
						stack_force.addBond(i, j, [0.36*unit.nanometer, 11.0*unit.kilocalorie_per_mole])
					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['ADE', 'THY']):
						stack_force.addBond(i, j, [0.36*unit.nanometer, 8.6*unit.kilocalorie_per_mole])
					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['ADE', 'CYT']):
						stack_force.addBond(i, j, [0.36*unit.nanometer, 7.8*unit.kilocalorie_per_mole])
					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['ADE', 'GUA']):
						stack_force.addBond(i, j, [0.36*unit.nanometer, 10.3*unit.kilocalorie_per_mole])
					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['THY', 'CYT']):
						stack_force.addBond(i, j, [0.36*unit.nanometer, 6.9*unit.kilocalorie_per_mole])
					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['THY', 'GUA']):
						stack_force.addBond(i, j, [0.36*unit.nanometer, 9.2*unit.kilocalorie_per_mole])
					if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['CYT', 'GUA']):
						stack_force.addBond(i, j, [0.36*unit.nanometer, 8.3*unit.kilocalorie_per_mole])
				elif row_i['chain'] == row_j['chain'] and abs(i-j) == 4:
					excluded_force.addBond(i, j)
				elif Counter([row_i['name'],row_j['name']]) != Counter(['ADE','THY']) and Counter([row_i['name'],row_j['name']]) != Counter(['CYT','GUA']):
						excluded_force.addBond(i, j)
				else:
					if Counter([row_i['name'],row_j['name']]) == Counter(['ADE','THY']):
						hb_force.addBond(i, j, [0.6*unit.nanometer,2.7*unit.kilocalorie_per_mole])
					elif Counter([row_i['name'],row_j['name']]) == Counter(['CYT','GUA']):
						hb_force.addBond(i, j, [0.55*unit.nanometer,3.3*unit.kilocalorie_per_mole])
					#else:
					#	excluded_force.addBond(i, j)
						
	stack_force.setUsesPeriodicBoundaryConditions(use_pbc)
	stack_force.setForceGroup(3)

	hb_force.setUsesPeriodicBoundaryConditions(use_pbc)
	hb_force.setForceGroup(4)

	excluded_force.setUsesPeriodicBoundaryConditions(use_pbc)
	excluded_force.setForceGroup(5)
	
	return stack_force, hb_force, excluded_force, vdw_force




def stack_12_10(dataframe, Mtype='DNA',use_pbc=True, force_group=3):
	"""
		U_stack = step(r_stack -r) * 
			\\delta_stack* [5 * (\\frac{\\sigma_stack}{r})^12 - 6 * (\\frac{\\sigma_{stack}}{r})^10]
	"""
	stack_expression = 'step(rstack -r) * (delta * (5* (sigma/r)^12 - 6*(sigma/r)^10))'
	stack_force      = mm.CustomBondForce(stack_expression)
	stack_force.addGlobalParameter('rstack', 0.62*unit.nanometer)

	stack_force.addPerBondParameter('sigma')
	stack_force.addPerBondParameter('delta')
	
	if Mtype == 'DNA':
		for i, row_i in dataframe.iterrows():
			for j, row_j in dataframe.iterrows():
				if i > j  and row_i['chain'] == row_j['chain']: # make sure they are in the same chain
					if abs(i-j) == 2:                           # make sure they are adjacent to each base
						if row_i['name'] == 'ADE' and row_j['name'] == 'ADE':
							stack_force.addBond(i, j, [0.36*unit.nanometer, 9.7*unit.kilocalorie_per_mole])
						if row_i['name'] == 'THY' and row_j['name'] == 'THY':
							stack_force.addBond(i, j, [0.36*unit.nanometer, 7.7*unit.kilocalorie_per_mole])
						if row_i['name'] == 'CYT' and row_j['name'] == 'CYT':
							stack_force.addBond(i, j, [0.36*unit.nanometer, 6.2*unit.kilocalorie_per_mole])
						if row_i['name'] == 'GUA' and row_j['name'] == 'GUA':
							stack_force.addBond(i, j, [0.36*unit.nanometer, 11.0*unit.kilocalorie_per_mole])
						if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['ADE', 'THY']):
							stack_force.addBond(i, j, [0.36*unit.nanometer, 8.6*unit.kilocalorie_per_mole])
						if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['ADE', 'CYT']):
							stack_force.addBond(i, j, [0.36*unit.nanometer, 7.8*unit.kilocalorie_per_mole])
						if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['ADE', 'GUA']):
							stack_force.addBond(i, j, [0.36*unit.nanometer, 10.3*unit.kilocalorie_per_mole])
						if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['THY', 'CYT']):
							stack_force.addBond(i, j, [0.36*unit.nanometer, 6.9*unit.kilocalorie_per_mole])
						if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['THY', 'GUA']):
							stack_force.addBond(i, j, [0.36*unit.nanometer, 9.2*unit.kilocalorie_per_mole])
						if row_i['name'] != row_j['name'] and Counter([row_i['name'],row_j['name']]) == Counter(['CYT', 'GUA']):
							stack_force.addBond(i, j, [0.36*unit.nanometer, 8.3*unit.kilocalorie_per_mole])
						
						
	stack_force.setUsesPeriodicBoundaryConditions(use_pbc)
	stack_force.setForceGroup(force_group)
	return stack_force


def HB_12_10(dataframe, Mtype='DNA', use_pbc=True, force_group=4):
	"""
		U_HB = step(r_HB -r) *
			\\delta_HB * [5* (\\frac{\\delta_HB}{r})^12 - 6* (\\frac{\\sigma_HB}{r})^10]
	"""
	hb_expression = """
				step(r_hb - r)*(delta_hb*(5*(sigma_hb/r)^12 - 6*(sigma_hb/r)^10))
	"""
	hb_force = mm.CustomBondForce(hb_expression)
	hb_force.addGlobalParameter('r_hb', 0.95*unit.nanometer)

	hb_force.addPerBondParameter('sigma_hb')
	hb_force.addPerBondParameter('delta_hb')

	if Mtype == 'DNA':
		for i, row_i in dataframe.iterrows():
			for j, row_j in dataframe.iterrows():
				if i > j and row_i['name'] != 'PO4' and row_j['name'] != 'PO4':
					if abs(i-j) > 4:
						if Counter([row_i['name'],row_j['name']]) == Counter(['ADE','THY']):
							hb_force.addBond(i, j, [0.6*unit.nanometer,2.7*unit.kilocalorie_per_mole])
						elif Counter([row_i['name'],row_j['name']]) == Counter(['CYT','GUA']):
							hb_force.addBond(i, j, [0.55*unit.nanometer,3.3*unit.kilocalorie_per_mole])

	hb_force.setUsesPeriodicBoundaryConditions(use_pbc)
	hb_force.setForceGroup(force_group)
	
	return hb_force

def excluded_12_6(dataframe, Mtype='DNA', use_pbc=True, force_group=5):
	"""
		U_excluded = step(r_exclude - r) * 4 * epsilon_ex * 
					[(sigma_ex/r)^12 - (sigma_ex/r)^6] + epsilon_ex
	"""
	excluded_expression = """
			step(rcut-r)*(4*epsilon*((sigma/r)^12 - (sigma/r)^6) + epsilon); rcut=2^(1/6)*sigma
	"""
	excluded_force      = mm.CustomBondForce(excluded_expression)
	excluded_force.addGlobalParameter('epsilon', 4*unit.kilocalorie_per_mole)
	excluded_force.addGlobalParameter('sigma', 0.55*unit.nanometer)

	# store all bonds info
	bonds = []
	for x, row_x in dataframe.iterrows():
		if len(row_x['bondlist']):
			for nxt in row_x['bondlist']:
				bonds.append([x, nxt])

	if Mtype == 'DNA':
		for i, row_i in dataframe.iterrows():
			for j, row_j in dataframe.iterrows():
				if i > j and row_i['chain'] == row_j['chain'] and abs(i-j) == 4:
					if row_i['name'] != 'PO4' and row_j['name'] != 'PO4':
						excluded_force.addBond(i, j)
				
				elif i > j and not any(set([i,j]) == set(bond) for bond in bonds):
					if Counter([row_i['name'],row_j['name']]) == Counter(['PO4','ADE']):
						excluded_force.addBond(i, j)
					elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','GUA']):
						excluded_force.addBond(i, j)
					elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','CYT']):
						excluded_force.addBond(i, j)
					elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','THY']):
						excluded_force.addBond(i, j)
					elif Counter([row_i['name'],row_j['name']]) == Counter(['PO4','PO4']):
						excluded_force.addBond(i, j)
						

	excluded_force.setUsesPeriodicBoundaryConditions(use_pbc)
	excluded_force.setForceGroup(force_group)
	
	return excluded_force

def vdw_12_6(dataframe, Mtype='DNA', use_pbc=True, force_group=6):
	"""
		U_vdw = step(r_cut - r) * (4 * epsilon * [(\\frac{\\sigma}{r})^12 - (\\frac{\\sigma}{r})^6] + epsilon* (1-\\lambda)) + \
				step(r - r_cut) * (4 * epsilon * [(\\frac{\\sigma}{r})^12 - (\\frac{\\sigma}{r})^6])
	"""

	vdw_expression = """
		step(r_vdw_cut - r)* (lj + (1-lambda)*epsilon) + step(r-r_vdw_cut)*(lambda*lj); \\
											lj = 4*epsilon*((sigma/r)^12 - (sigma/r)^6); \\
											r_vdw_cut = 2^(1/6)*sigma; \\
											lambda = 0.5*(lambda1 + lambda2); \\
											sigma=0.5*(sigma1+sigma2)
	"""

	vdw_force           = mm.CustomNonbondedForce(vdw_expression)
	vdw_force.addPerParticleParameter('epsilon')
	vdw_force.addPerParticleParameter('sigma')
	vdw_force.addPerParticleParameter('lambda')

	# store all bonds info
	bonds = []
	for x, row_x in dataframe.iterrows():
		if len(row_x['bondlist']):
			for nxt in row_x['bondlist']:
				bonds.append((x, nxt))

	# assign parameters
	max_radii = 0
	for i, row in dataframe.iterrows():
		if Mtype == 'DNA':
			eps   = 0.2
			radii = HPS_DNA.loc[HPS_DNA['atomtype']==row['name'],'sigma'].values[0]
			lba = HPS_DNA.loc[HPS_DNA['atomtype']==row['name'],'lambda'].values[0]
		elif Mtype == 'protein':
			eps  = 0.2
			radii = HPS_PROT.loc[HPS_PROT['atomtype']==row['ATYPE'],'sigma'].values[0]
			lba   = HPS_PROT.loc[HPS_PROT['atomtype']==row['ATYPE'],'lambda'].values[0]
		if radii > max_radii:
			max_radii = radii
		# three parameters are required:
		# epsilon
		# sigma
		# lambda
		vdw_force.addParticle([esp* unit.kilocalorie_per_mole, radii* unit.nanometer, lba])
		
	
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

	if Mtype == 'DNA':
		for i, row in dataframe.iterrows():
			chg = HPS_DNA.loc[HPS_DNA['atomtype']==row['name'],'charge'].values[0]
			elec_force.addParticle([chg])
		#	if row['name'] == 'PO4':
		#		elec_force.addParticle([-1])
		#	else:
		#		elec_force.addParticle([0])

	elif Mtype == 'protein':
		for i, row in dataframe.iterrows():
			chg = HPS_PROT.loc[HPS_PROT['atomtype']==row['name'],'charge'].values[0]
			elec_force.addParticle([chg])
			
	# store all bonds info
	bonds = []
	for x, row_x in dataframe.iterrows():
		if len(row_x['bondlist']):
			for nxt in row_x['bondlist']:
				bonds.append((x, nxt))

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

