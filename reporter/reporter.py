#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys, os
from pathlib import Path
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import time

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

	# add forces/function forms'
	print('Starting add forces: Bond')
	s_time = time.perf_counter()
	U_bond = harmonic_bond_form(data, Mtype=Mtype)
	e_time = time.perf_counter()
	print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	if 'DNA' in Mtypes:
		print('Starting add forces: Angle')
		s_time = time.perf_counter()
		U_angle = harmonic_cos_angle_form(data, Mtype=Mtype)
		e_time = time.perf_counter()
		print('Elasping: {:.3f} seconds'.format(e_time - s_time))		

		print('Starting add forces: stack + HB + Exclude')
		s_time = time.perf_counter()
		U_stack, U_HB, U_exclude, U_vdw_DNA = stack_hb_exclude_DNA(data, Mtype=Mtype)
		e_time = time.perf_counter()
		print('Elasping: {:.3f} seconds'.format(e_time - s_time))
	#elif 'protein' in Mtypes:
	
	if 'protein' in Mtypes:
		print('Starting add forces: vdW')
		s_time = time.perf_counter()
		U_vdw  = vdw_12_6(data, Mtype=Mtype)
		e_time = time.perf_counter()
		print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	print('Starting add forces: Coul')
	s_time = time.perf_counter()
	U_elec = elec(data, Mtype=Mtype, dielectric=conf['dielectric'])
	e_time = time.perf_counter()
	print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	# using defined forces
	system.addForce(U_bond)
	if 'DNA' in Mtypes:
		system.addForce(U_angle)
		system.addForce(U_stack)
		system.addForce(U_HB)
		system.addForce(U_exclude)
		system.addForce(U_vdw_DNA)

	if 'protein' in Mtypes:
		system.addForce(U_vdw)
	system.addForce(U_elec)

	TEMP = conf['temp']
	plat = mm.Platform.getPlatformByName(conf['platform'])
	prop = conf['prop']
	ts   = conf['timestep']

	integrator = mm.LangevinMiddleIntegrator(TEMP*unit.kelvin, 0.01/unit.picosecond, ts*unit.picoseconds)
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

def NVT(data, conf):

	print(data)
	topology, system, = build_system(data)
	print('Done with system construction!')
	
	Mtypes = ['protein']*conf['protein'] + ['DNA']*conf['DNA']
	print(Mtypes)


	vectors = (conf['pbc']) * unit.nanometers
	system.setDefaultPeriodicBoxVectors(*vectors)
	topology.setPeriodicBoxVectors(vectors)

	# load the coordinate you want to estimate the energy
	# convert series into numpy for performance
	xyz = np.stack(data['xyz'].values) * 0.1 * unit.nanometer
	
	# add forces/function forms
	print('Starting add forces: Bond')
	s_time = time.perf_counter()
	U_bond = harmonic_bond_form(data)
	e_time = time.perf_counter()
	print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	if 'DNA' in Mtypes:
		print('Starting add forces: Angle')
		s_time = time.perf_counter()
		U_angle = harmonic_cos_angle_form(data)
		e_time = time.perf_counter()
		print('Elasping: {:.3f} seconds'.format(e_time - s_time))

		print('Starting add forces: Stack + HB + Exclude + vdW')
		s_time = time.perf_counter()
		U_stack, U_HB, U_exclude, U_vdw_DNA = stack_hb_exclude_DNA(data)
		e_time = time.perf_counter()
		print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	if 'protein' in Mtypes:
		print('Starting add forces: vdW')
		s_time = time.perf_counter()
		U_vdw = vdw_12_6(data)
		e_time = time.perf_counter()
		print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	print('Starting add forces: Coul')
	s_time = time.perf_counter()
	U_elec = elec(data, dielectric=conf['dielectric'])
	e_time = time.perf_counter()
	print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	# using defined forces
	system.addForce(U_bond)
	if 'DNA' in Mtypes:
		system.addForce(U_angle)
		system.addForce(U_stack)
		system.addForce(U_HB)
		system.addForce(U_exclude)
		system.addForce(U_vdw_DNA)

	if 'protein' in Mtypes:
		system.addForce(U_vdw)
	system.addForce(U_elec)

	TEMP = conf['temperature']
	plat = mm.Platform.getPlatformByName(conf['platform'])
	prop = conf['prop']
	ts   = conf['timestep']

	#print('DEBUG:',ts)

	integrator = mm.LangevinMiddleIntegrator(TEMP*unit.kelvin, 0.01/unit.picosecond, ts*unit.picoseconds)
	simulation = app.Simulation(topology, system, integrator, plat, prop)
	simulation.context.setPositions(xyz)
	
	# section for minimization
	print('Starting monimization !')
	s_time = time.perf_counter()
	simulation.minimizeEnergy()
	positions = simulation.context.getState(getPositions=True).getPositions()
	with open('minimized_structure.pdb', 'w') as f:
		app.PDBFile.writeFile(topology, positions, f)	
	e_time = time.perf_counter()
	print('Elasping: {:.3f} seconds'.format(e_time - s_time))

	simulation.context.setVelocitiesToTemperature(TEMP * unit.kelvin)
	simulation.reporters.append(app.DCDReporter(conf['trajfile'], conf['interval'], enforcePeriodicBox=False)) # enforcePeriodicBox
																												# True = wrap; False = leave it there
	
	labels = {1:'bond',2:'angle',3:'stacking',4:'hb',5:'excluding',6:'vdw',7:'elec'}
	
	energies = {}
	for gid, name in labels.items():
		state = simulation.context.getState(getEnergy=True, groups={gid})
		energies[name] = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

	# chkpoint reporter
	simulation.reporters.append(app.CheckpointReporter(conf['rstfile'],reportInterval=conf['interval']))

	# energy reporter
	simulation.reporters.append(EnergyDecompositionReporter(conf['log'],reportInterval=conf['interval'],labels=labels))
	simulation.step(conf['step'])
	simulation.saveCheckpoint(conf['rstfile'])
	return None


