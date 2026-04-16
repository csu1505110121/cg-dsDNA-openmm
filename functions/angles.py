#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os, sys
import openmm as mm
import openmm.app as app
import openmm.unit as unit

def harmonic_cos_angle_form(dataframe, Mtype='DNA', use_pbc=True, force_group=2):
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

	#### old script
##	for i, row in dataframe.iterrows():
##		Mtype = row['Mtype']
##		if Mtype == 'DNA':
##			atomidx = row['anglelist']
##			if len(atomidx):
##				angle_force.addAngle(atomidx[0], i, atomidx[1], [180*unit.degree, 40*unit.kilocalorie_per_mole])

	angle_force.setUsesPeriodicBoundaryConditions(use_pbc)
	angle_force.setForceGroup(force_group)
	return angle_force



