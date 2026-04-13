#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys, os
from pathlib import Path
import openmm as mm
import openmm.app as app
import openmm.unit as unit

current_file = Path(__file__).resolve()
project_root = current_file.parent.parent

if str(project_root) not in sys.path:
	sys.path.append(str(project_root))

from tools.model_build import build_system
from tools.converting import readDUMPFILE
from parameters.load_params import readParam
from functions.angles import harmonic_cos_angle_form
from functions.bonds import harmonic_bond_form
from functions.nonbonded import *


class EnergyDecompositionReporter(object):
	def __init__(self, file, reportInterval, labels):
		self._out = open(file, 'w')
		self._interval = reportInterval
		self._labels = labels
		self._out.write("Step " + " ".join(labels.values()) + "\n")

	def describeNextReport(self, simulation):
		steps_until_next_report = self._interval - simulation.currentStep % self._interval
		return (steps_until_next_report, False, False, False, True)

	def report(self, simulation, state):
		step = simulation.currentStep
		self._out.write(str(step))

		for g in self._labels:
			s = simulation.context.getState(
				getEnergy=True,
				groups={g}
			)
			e = s.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
			self._out.write(f" {e}")

		self._out.write("\n")
		self._out.flush()


def _sp_estimation(data, conf,  Mtype='DNA'):
	"""
		single point estimator
	"""
	topology, system, = build_system(data, Mtype=Mtype)

	Mtypes = ['protein']*conf['protein'] + ['DNA']*conf['DNA']

	vectors = (conf['pbc']) * unit.nanometers 
	
	system.setDefaultPeriodicBoxVectors(*vectors)
	topology.setPeriodicBoxVectors(vectors)

	# load the coordinate you want to estimate the energy
	xyz = readDUMPFILE(conf['xyz']) * 0.1 * unit.nanometer

	# add forces/function forms
	U_bond = harmonic_bond_form(data, Mtype=Mtype)
	if 'DNA' in Mtypes:
		U_angle = harmonic_cos_angle_form(data, Mtype=Mtype)
		U_stack, U_HB, U_exclude, U_vdw_DNA = stack_hb_exclude_DNA(data, Mtype=Mtype)
	#elif 'protein' in Mtypes:
	
	U_vdw  = vdw_12_6(data, Mtype=Mtype)
	U_elec = elec(data, Mtype=Mtype, dielectric=conf['dielectric'])

	# using defined forces
	system.addForce(U_bond)
	if 'DNA' in Mtypes:
		system.addForce(U_angle)
		system.addForce(U_stack)
		system.addForce(U_HB)
		system.addForce(U_exclude)
		system.addForce(U_vdw_DNA)

	system.addForce(U_vdw)
	system.addForce(U_elec)

	TEMP = conf['temp']
	plat = mm.Platform.getPlatformByName(conf['platform'])
	prop = conf['prop']

	integrator = mm.LangevinMiddleIntegrator(TEMP*unit.kelvin, 0.01/unit.picosecond, 0.010*unit.picoseconds)
	simulation = app.Simulation(topology, system, integrator, plat, prop)
	simulation.context.setPositions(xyz)

	simulation.context.setVelocitiesToTemperature(TEMP * unit.kelvin)

	labels = {1:'bond',2:'angle',3:'stacking',4:'hb',5:'excluding',6:'vdw',7:'elec'}
	
	energies = {}
	for gid, name in labels.items():
		state = simulation.context.getState(getEnergy=True, groups={gid})
		energies[name] = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

	simulation.step(conf['step'])

	return energies

def NVT(data, conf, Mtype='DNA'):
	topology, system, = build_system(data, Mtype=Mtype)
	
	vectors = (conf['pbc']) * unit.nanometers
	system.setDefaultPeriodicBoxVectors(*vectors)
	topology.setPeriodicBoxVectors(vectors)

	# load the coordinate you want to estimate the energy
	xyz = data['xyz'].values * 0.1 * unit.nanometer
	
	# add forces/function forms
	U_bond = harmonic_bond_form(data, Mtype=Mtype)
	U_angle = harmonic_cos_angle_form(data, Mtype=Mtype)
	U_stack, U_HB, U_exclude, U_vdw = stack_hb_exclude_DNA(data, Mtype=Mtype)
	U_elec = elec(data, Mtype=Mtype, dielectric=conf['dielectric'])

	# using defined forces
	system.addForce(U_bond)
	system.addForce(U_angle)
	system.addForce(U_stack)
	system.addForce(U_HB)
	system.addForce(U_exclude)
	system.addForce(U_vdw)
	system.addForce(U_elec)

	TEMP = conf['temp']
	plat = mm.Platform.getPlatformByName(conf['platform'])
	prop = conf['prop']

	integrator = mm.LangevinMiddleIntegrator(TEMP*unit.kelvin, 0.01/unit.picosecond, 0.010*unit.picoseconds)
	simulation = app.Simulation(topology, system, integrator, plat, prop)
	simulation.context.setPositions(xyz)
	
	simulation.minimizeEnergy()
	
	simulation.context.setVelocitiesToTemperature(TEMP * unit.kelvin)
	simulation.reporters.append(app.DCDReporter('traj.dcd', 10_000, enforcePeriodicBox=False))
	
	labels = {1:'bond',2:'angle',3:'stacking',4:'hb',5:'excluding',6:'vdw',7:'elec'}
	
	energies = {}
	for gid, name in labels.items():
		state = simulation.context.getState(getEnergy=True, groups={gid})
		energies[name] = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

	simulation.reporters.append(EnergyDecompositionReporter(conf['log'],reportInterval=10_000,labels=labels))
	simulation.step(conf['step'])
	return None


