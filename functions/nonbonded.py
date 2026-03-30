#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd


def stack_12_10(use_pbc, force_group=3):
	"""
		U_stack = step(r_stack -r) * 
			\delta_stack* [5 * (\frac{\sigma_stack}{r})^12 - 6 * (\frac{\sigma_{stack}}{r})^10]
	"""
	stack_expression = 'step(rstack -r) * (delta * (5* (sigma/r)^12 - 6*(sigma/r)^10))'
	stack_force      = mm.CustomBondForce(stack_expression)
	stack_force.addGlobalParameter('rstack', 0.62*unit.nanometer)

	stack_force.addPerBondParameter('sigma')
	stack_force.addPerBondParameter('delta')

	stack_force.setUsesPeriodicBoundaryConditions(use_pbc)
	stack_force.setForceGroup(force_group)
	return stack_force


def HB_12_10(use_pbc, force_group=4):
	"""
		U_HB = step(r_HB -r) *
			\delta_HB * [5* (\frac{\delta_HB}{r})^12 - 6* (\frac{\sigma_HB}{r})^10]
	"""
	hb_expression = """
				step(r_hb - r)*(delta_hb*(5*(sigma_hb/r)^12 - 6*(sigma_hb/r)^10))
	"""
	hb_force = mm.CustomBondForce(hb_expression)
	hb_force.addGlobalParameter('r_hb', 0.95*unit.nanometer)

	hb_force.addPerBondParameter('sigma_hb')
	hb_force.addPerBondParameter('delta_hb')

	hb_force.setUsesPeriodicBoundaryConditions(use_pbc)
	hb_force.setForceGroup(force_group)
	
	return hb_force

def excluded_12_6(use_pbc, force_group=5):
	"""
		U_excluded = step(r_exclude - r) * 4 * epsilon_ex * 
					[(sigma_ex/r)^12 - (sigma_ex/r)^6] + epsilon_ex
	"""
	excluded_expression = """
			step(rcut-r)*(4*epsilon*((sigma/r)^12 - (sigma/r)^6) + epsilon); rcut=2^(1/6)*sigma
	"""
	excluded_force      = mm.CustomBondedForce(excluded_expression)
	excluded_force.addGlobalParameter('epsilon', 4*unit.kilocalorie_per_mole)
	excluded_force.addGlobalParameter('sigma', 0.55*unit.nanometer)

	excluded_force.setUsesPeriodicBoundaryConditions(use_pbc)
	excluded_force.setForceGroup(force_group)
	
	return excluded_force

def vdw_12_6(use_pbc, force_group=6):
	"""
		U_vdw = step(r_cut - r) * (4 * epsilon * [(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6] + epsilon* (1-\lambda)) + \
				step(r - r_cut) * (4 * epsilon * [(\frac{\sigma}{r})^12 - (\frac{\sigma}{r})^6])
	"""

	vdw_expression = """
		step(r_vdw_cut - r)*step(r_vdw-r) * (4*epsilon*((sigma/r)^12 - (sigma/r)^6) + (1-lambda)*epsilon)
	"""
	vdw_force           = mm.CustomBondedForce(vdw_expression)
	vdw_force.addGlobalParameter('lambda',0)
	vdw_force.addGlobalParameter('r_vdw',0.55*2**(1/6)*unit.nanometer)
	vdw_force.addGlobalParameter('r_vdw_cut',0.6174*unit.nanometer)

	vdw_force.addGlobalParameter('epsilon',4*unit.kilocalorie_per_mole)
	vdw_force.addGlobalParameter('sigma',0.55*unit.nanometer)

	vdw_force.setUsesPeriodicBoundaryConditions(use_pbc)
	vdw_force.setForceGroup(force_group)

	return vdw_force


def elec(use_pbc, force_group=7):
	"""
		U_electrostatics = \frac{q_i * q_j}{4*\pi*\epsilon_0 * \epsilon_r * r_{ij}} \exp^{-r_ij/\kappa}
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
	elec_force.addGlobalParameter('dielectric_water',80)

	elec_force.addPerParticleParameter('q')

	for atom in topology.atoms():
		pass

	elec_force.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()],1)

	elec_force.setUsesPeriodicBoundaryConditions(use_pbc)
	elec_force.setCutoffDistance(3.5*unit.nanometer)
	elec_force.setForceGroup(force_group)

	return elec_force





