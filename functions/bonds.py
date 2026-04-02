#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os,sys
import openmm as mm
import openmm.app as app
import openmm.unit as unit

def harmonic_bond_form(dataframe, Mtype='protein', use_pbc=True, force_group=1):
	"""
	harmonic bond potential
	U_bond = \sum K_{bond} (r- r_0)^2

	K_bond -> K: force constant in the unit of kcal/mol/A2
	r_0    -> r0: equilibrium bond length in the unit of nm

	-Mtype: protein | DNA


	"""
	bond_expression = 'K*(r - r0)^2'
	bond_force      = mm.CustomBondForce(bond_expression)
	bond_force.addPerBondParameter('r0')
	bond_force.addPerBondParameter('K')

	if Mtype == 'protein':
		for i, row in dataframe.iterrows():
			atomidx = row['bondlist'][0]
			bond_force.addBond(i, atomidx, [0.38*unit.nanometer, 10*unit.kilocalorie_per_mole/unit.angstrom**2])

	elif Mtype == 'DNA':
		for i, row in dataframe.iterrows():
			atomidx = row['bondlist'][]
			if len(atomidx) == 2:
				# the first index is always for BASE connecting with PO4
				bond_force.addBond(i, atomidx[0], [0.575*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
				# the second index is for PO4
				bond_force.addBond(i, atomidx[1], [0.55*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
			elif len(atomidx) == 1:
				bond_force.addBond(i, atomidx[0], [0.575*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
	elif Mtype == 'organic':
		pass

	bond_force.setUsesPeriodicBoundaryConditions(use_pbc)
	bond_force.setForceGroup(force_group)
	return bond_force


