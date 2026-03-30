#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os,sys


def harmonic_bond_form(use_pbc, force_group=1):
	"""
	harmonic bond potential
	U_bond = \sum K_{bond} (r- r_0)^2

	K_bond -> K: force constant in the unit of kcal/mol/A2
	r_0    -> r0: equilibrium bond length in the unit of nm

	"""
	bond_expression = 'K*(r - r0)^2'
	bond_force      = mm.CustomBondForce(bond_expression)
	bond_force.addPerBondParameter('r0')
	bond_force.addPerBondParameter('K')

	for bond in topology.bonds():
		bond_force.addBond(bond.atom1.index, bond.atom2.index, [])

	bond_force.setUsesPeriodicBoundaryConditions(use_pbc)
	bond_force.setForceGroup(force_group)
	return bond_force
