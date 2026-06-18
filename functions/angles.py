#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os, sys
import openmm as mm
import openmm.app as app
import openmm.unit as unit

def harmonic_cos_angle_form01(dataframe, Mtype='DNA', use_pbc=True, force_group=2):
	"""
	cosine/squared-style harmonic potential:
	U_angle = \\sum K_angle (\\cos \\theta - \\cos \\theta_0)^2

	K_angle: force constant in the unit of kcal/mol
	theta0: equilibrium angle in the unit of rad
	"""
	
	angle_expression = 'K * (cos(theta) - cos(theta0))^2'
	angle_force      = mm.CustomAngleForce(angle_expression)
	angle_force.addPerAngleParameter("theta0")
	angle_force.addPerAngleParameter("K")


	# only need to handle DNA when Mtype == 'DNA'

	dna_df = dataframe[dataframe['Mtype']=='DNA']
	dna_indices = dna_df.index.tolist()
	dna_anglelist = dna_df['anglelist'].tolist()


	for _indice, _anglelist in zip(dna_indices, dna_anglelist):
		if len(_anglelist):
			angle_force.addAngle(_anglelist[0], _indice, _anglelist[1], [180*unit.degree, 40*unit.kilocalorie_per_mole])


	angle_force.setUsesPeriodicBoundaryConditions(use_pbc)
	angle_force.setForceGroup(force_group)
	return angle_force





def harmonic_cos_angle_form(dataframe, use_pbc=True, force_group=2):
	"""
	cosine/squared-style harmonic potential:
	U_angle = \\sum K_angle (\\cos \\theta - \\cos \\theta_0)^2

	K_angle: force constant in the unit of kcal/mol
	theta0: equilibrium angle in the unit of rad
	"""
	
	angle_expression = 'K * (cos(theta) - cos(theta0))^2'
	angle_force      = mm.CustomAngleForce(angle_expression)
	angle_force.addPerAngleParameter("theta0")
	angle_force.addPerAngleParameter("K")


	_names = dataframe['name'].to_dict()
	_Mtypes= dataframe['Mtype'].tolist()
	_atypes= dataframe['ATYPE'].tolist()
	_anglelists=dataframe['anglelist'].tolist()
	_indices  =dataframe.index.tolist()

	#print(_Mtypes, _atypes, _names, _anglelists)


	for _indice, _mtype, _atype, _anglelist in zip(_indices, _Mtypes, _atypes, _anglelists):
		if len(_anglelist):
			if _mtype == 'DNA':
				angle_force.addAngle(_anglelist[0], _indice, _anglelist[1], [180*unit.degree, 40*unit.kilocalorie_per_mole])
			elif _mtype == 'nucleotide':
				if _atype in ['ATP', 'CTP', 'GTP', 'TTP']:
					if _names[_indice] == 'P1' and _names[_anglelist[0]] == 'P2' and _names[_anglelist[1]] == 'P1':  # angle for P1-P1-P2
						angle_force.addAngle(_anglelist[0], _indice, _anglelist[1], [120*unit.degree, 40*unit.kilocalorie_per_mole])
					elif _names[_indice] == 'P1' and _names[_anglelist[0]] == 'P1' and _names[_anglelist[1]] == 'CR':
						angle_force.addAngle(_anglelist[0], _indice, _anglelist[1], [90*unit.degree, 40*unit.kilocalorie_per_mole])
					elif _names[_indice] == 'CR' and _names[_anglelist[0]] == 'P1' and _names[_anglelist[1]] in ['XA', 'XG', 'XC', 'XT']:
						angle_force.addAngle(_anglelist[0], _indice, _anglelist[1], [90*unit.degree, 40*unit.kilocalorie_per_mole])


	angle_force.setUsesPeriodicBoundaryConditions(use_pbc)
	angle_force.setForceGroup(force_group)
	return angle_force