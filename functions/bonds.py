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

	_names = dataframe['name'].to_dict()
	_Mtypes= dataframe['Mtype'].tolist()
	_atypes= dataframe['ATYPE'].tolist()
	_bondlists=dataframe['bondlist'].tolist()
	_indices  =dataframe.index.tolist()

	for _indice, _mtype, _bondlist in zip(_indices, _Mtypes, _bondlists):
		if len(_bondlist):
			if _mtype == 'protein':
				for nxt in _bondlist:
					bond_force.addBond(_indice, nxt, [0.38*unit.nanometer, 10*unit.kilocalorie_per_mole/unit.angstrom**2])
			elif _mtype == 'DNA':
				for nxt in _bondlist:
					_name = _names[nxt]
					if _name == 'PO4': # bond for PO4-PO4
						bond_force.addBond(_indice, nxt, [0.55*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
					elif _name != 'PO4': # bond for PO4-BASE
						bond_force.addBond(_indice, nxt, [0.575*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
			elif _mtype == 'nucleotide':
				for nxt in _bondlist:
					_name = _names[nxt]
					if _name == 'P1' and _names[_indice] == 'P2': # bond for PG-PB
						bond_force.addBond(_indice, nxt, [0.35*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
					elif _name == 'P1' and _names[_indice] == 'P1': # bond for PB-PA
						bond_force.addBond(_indice, nxt, [0.35*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
					elif _name == 'CR' and _names[_indice] == 'P1': # bond for PA-CR
						bond_force.addBond(_indice, nxt, [0.45*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
					elif _name in ['XA', 'XG', 'XC', 'XT'] and _names[_indice] == 'CR': # bond for PO4-BASE
						bond_force.addBond(_indice, nxt, [0.435*unit.nanometer, 50*unit.kilocalorie_per_mole/unit.angstrom**2])
			else:
				raise('In applying Bond force, however the molecule type {} is not implemented yet'.format(_mtype))
	


	bond_force.setUsesPeriodicBoundaryConditions(use_pbc)
	bond_force.setForceGroup(force_group)
	return bond_force


