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

	#if Mtype == 'protein':
	for i, row in dataframe.iterrows():
		Mtype = row['Mtype']
		if len(row['bondlist']):
			if Mtype == 'protein':
				atomidx = row['bondlist'][0]
				bond_force.addBond(i, atomidx, [0.38*unit.nanometer, 10*unit.kilocalorie_per_mole/unit.angstrom**2])
			elif Mtype == 'DNA':
				atomidx = row['bondlist']
				if len(atomidx):
					for nxt in atomidx:
						if dataframe.iloc[nxt]['name'] == 'PO4':
							bond_force.addBond(i, nxt, [0.55*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
						elif dataframe.iloc[nxt]['name'] != 'PO4':
							bond_force.addBond(i, nxt, [0.575*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
			elif Mtype == 'organic':
				pass
			else:
				raise('In applying Bond force, however the molecule type {} is not implemented yet'.format(Mtype))

#	elif Mtype == 'DNA':
#		for i, row in dataframe.iterrows():
#			atomidx = row['bondlist']
#			if len(atomidx):
#				for nxt in atomidx:
#					if dataframe.iloc[nxt]['name'] == 'PO4':
#						bond_force.addBond(i, nxt, [0.55*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
#					elif dataframe.iloc[nxt]['name'] != 'PO4':
#						bond_force.addBond(i, nxt, [0.575*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
#	elif Mtype == 'organic':
#		pass

	bond_force.setUsesPeriodicBoundaryConditions(use_pbc)
	bond_force.setForceGroup(force_group)
	return bond_force


