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
from tools.addRigid import _addrigid
from tools.io import STRSUMMARY
from parameters.load_params import readParam
from functions.angles import harmonic_cos_angle_form
from functions.bonds import harmonic_bond_form
from functions.nonbonded import *

class BoxSizeReporter(object):
    def __init__(self, file, reportInterval):
        # file 可以是文件名字符串，也可以是 sys.stdout (控制台)
        self._file = open(file, 'w') if isinstance(file, str) else file
        self._reportInterval = reportInterval

    def __del__(self):
        if hasattr(self, '_file') and isinstance(self._file, str):
            self._file.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, True) # 最后一个 True 获取 PBC 向量

    def report(self, simulation, state):
        # 获取当前状态的周期性盒子向量
        # vectors 是一个包含三个 Vec3 对象的元组
        vectors = state.getPeriodicBoxVectors()
        
        # 将单位转换为纳米并提取数值
        v1 = vectors[0].value_in_unit(unit.nanometer)
        v2 = vectors[1].value_in_unit(unit.nanometer)
        v3 = vectors[2].value_in_unit(unit.nanometer)
        
        step = simulation.currentStep
        
        # 格式化打印：打印出 3x3 矩阵的形式或单行形式
        output = (
            f"Step: {step}\n"
            f"  v1: ({v1[0]:10.4f}, {v1[1]:10.4f}, {v1[2]:10.4f}) nm\n"
            f"  v2: ({v2[0]:10.4f}, {v2[1]:10.4f}, {v2[2]:10.4f}) nm\n"
            f"  v3: ({v3[0]:10.4f}, {v3[1]:10.4f}, {v3[2]:10.4f}) nm"
        )
        
        print(output, file=self._file)
        # 如果是打印到文件，建议 flush 保证实时写入
        if hasattr(self._file, 'flush'):
            self._file.flush()


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


def _sp_estimation(data, conf):
	"""
		single point estimator
	"""
	topology, system, = build_system(data)


	Mtypes = ['protein']*conf['protein'] + ['DNA']*conf['DNA']

	vectors = (conf['pbc']) * unit.nanometers 
	
	system.setDefaultPeriodicBoxVectors(*vectors)
	topology.setPeriodicBoxVectors(vectors)

	# load the coordinate you want to estimate the energy
	xyz = readDUMPFILE(conf['xyz']) * 0.1 * unit.nanometer

	# add forces/function forms'
	print('Starting add forces: Bond')
	s_time = time.perf_counter()
	U_bond = harmonic_bond_form(data)
	e_time = time.perf_counter()
	print('Elasped: {:.3f} seconds'.format(e_time - s_time))

	if 'DNA' in Mtypes:
		print('Starting add forces: Angle')
		s_time = time.perf_counter()
		U_angle = harmonic_cos_angle_form(data)
		e_time = time.perf_counter()
		print('Elasped: {:.3f} seconds'.format(e_time - s_time))		

		print('Starting add forces: stack + HB + Exclude')
		s_time = time.perf_counter()
		U_stack, U_HB, U_exclude, U_vdw_DNA = stack_hb_exclude_DNA(data)
		e_time = time.perf_counter()
		print('Elasped: {:.3f} seconds'.format(e_time - s_time))
	#elif 'protein' in Mtypes:
	
	if 'protein' in Mtypes:
		print('Starting add forces: vdW')
		s_time = time.perf_counter()
		U_vdw  = vdw_12_6(data)
		e_time = time.perf_counter()
		print('Elasped: {:.3f} seconds'.format(e_time - s_time))

	print('Starting add forces: Coul')
	s_time = time.perf_counter()
	U_elec = elec(data, dielectric=conf['dielectric'], C_salt=conf['C_salt'])
	e_time = time.perf_counter()
	print('Elasped: {:.3f} seconds'.format(e_time - s_time))

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


def _get_xyz_from_checkpoint(topology, system, conf):
    """从checkpoint文件中读取坐标"""
    # 用一个临时的integrator和simulation来读取checkpoint
    rstfile_p = conf['rstfile_p']
    tmp_integrator = mm.VerletIntegrator(0.001)
    tmp_plat = mm.Platform.getPlatformByName(conf['platform'])  # 用CPU读，避免GPU的问题
    tmp_prop = conf['prop']
    tmp_sim = app.Simulation(topology, system, tmp_integrator, tmp_plat, tmp_prop)
    
    with open(rstfile_p, 'rb') as f:
        tmp_sim.loadCheckpoint(f)
    
    state = tmp_sim.context.getState(getPositions=True)
    positions = state.getPositions()
    
    del tmp_sim
    return positions


def md(data, conf):

	# store the all infomation 
	sys_info = STRSUMMARY(data)

	#print(data)
	topology, system, = build_system(data)
	print('Done with system construction!')
	
	Mtypes = ['protein']*conf['protein'] + ['DNA']*conf['DNA']
	#print(Mtypes)


	vectors = (conf['pbc']) * unit.nanometers
	system.setDefaultPeriodicBoxVectors(*vectors)
	topology.setPeriodicBoxVectors(vectors)

	# load the coordinate you want to estimate the energy
	# convert series into numpy for performance
	#if conf['md'] == 'md': # this is for initializing the simulation
	#	xyz = np.stack(data['xyz'].values) * 0.1 * unit.nanometer
	
	# add forces/function forms
	print('Starting add forces: Bond')
	s_time = time.perf_counter()
	U_bond = harmonic_bond_form(data)
	e_time = time.perf_counter()
	print('Elasped: {:.3f} seconds'.format(e_time - s_time))

	if 'DNA' in Mtypes:
		print('Starting add forces: Angle')
		s_time = time.perf_counter()
		U_angle = harmonic_cos_angle_form(data)
		e_time = time.perf_counter()
		print('Elasped: {:.3f} seconds'.format(e_time - s_time))

		print('Starting add forces: Stack + HB + Exclude + vdW')
		s_time = time.perf_counter()
		U_stack, U_HB, U_exclude, U_vdw_DNA = stack_hb_exclude_DNA(data)
		e_time = time.perf_counter()
		print('Elasped: {:.3f} seconds'.format(e_time - s_time))

	if 'protein' in Mtypes:
		print('Starting add forces: vdW')
		s_time = time.perf_counter()
		U_vdw = vdw_12_6(data)
		e_time = time.perf_counter()
		print('Elasped: {:.3f} seconds'.format(e_time - s_time))

	print('Starting add forces: Coul')
	s_time = time.perf_counter()
	U_elec = elec(data, dielectric=conf['dielectric'], C_salt=conf['C_salt'])
	e_time = time.perf_counter()
	print('Elasped: {:.3f} seconds'.format(e_time - s_time))

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


	# remove the CMMotion
	force  = mm.CMMotionRemover()
	system.addForce(force)

	TEMP = conf['temperature']
	plat = mm.Platform.getPlatformByName(conf['platform'])
	prop = conf['prop']
	ts   = conf['timestep']

	#simulation = app.Simulation(topology, system, integrator, plat, prop)

	if conf['setrigid'] != 0:
		if conf['md'] == 'md':
			xyz = np.stack(data['xyz'].values) * 0.1 * unit.nanometer
		else:
			xyz = _get_xyz_from_checkpoint(topology, system, conf)

		system = _addrigid(system, xyz, sys_info)


	if conf['ensemble'] == 'NPT':
		print(f"Adding MonteCarloBarostat at {conf['ref_bar']} bar")
		ref_bar = float(conf['ref_bar'])
		r_bar_freq = int(conf['ref_bar_freq'])
		system.addForce(mm.MonteCarloBarostat(ref_bar *unit.bar, TEMP * unit.kelvin, r_bar_freq))


	integrator = mm.LangevinMiddleIntegrator(TEMP*unit.kelvin, 0.01/unit.picosecond, ts*unit.picoseconds)

	simulation = app.Simulation(topology, system, integrator, plat, prop)

	if conf['md'] == 'md':
		xyz = np.stack(data['xyz'].values) * 0.1 * unit.nanometer
		simulation.context.setPositions(xyz)
		print('Starting minimization...')
		simulation.minimizeEnergy()

		positions = simulation.context.getState(getPositions=True).getPositions()
		with open('minimized_structure.pdb', 'w') as f:
			app.PDBFile.writeFile(topology, positions, f)	
		
		simulation.context.setVelocitiesToTemperature(TEMP * unit.kelvin)

	elif conf['md'] in ['restart', 'slab']:
		with open(conf['rstfile_p'],'rb') as f:
			simulation.loadCheckpoint(f)

		if conf['md'] == 'slab':
			state = simulation.context.getState(getPositions=True, getVelocities=True)
			old_vecs = state.getPeriodicBoxVectors()
			new_vecs = (old_vecs[0], old_vecs[1], old_vecs[2] * conf['vecScale'])
			simulation.context.setPeriodicBoxVectors(*new_vecs)
			simulation.minimizeEnergy()
























#	if conf['md'] == 'md':
#		#simulation.context.setPositions(xyz)
#		# set up rigid bodies
#		if conf['setrigid'] != 0:
#			#state = simulation.context.getState(getPositions=True)
#			#positions = state.getPositions()
#			print('Applying constraints!')
#			system = _addrigid(system, xyz, sys_info)
#			print("After _addrigid vectors:", system.getDefaultPeriodicBoxVectors())
#
#
#
#
#			system.setDefaultPeriodicBoxVectors(*vectors)
#			print("After setDefault vectors:", system.getDefaultPeriodicBoxVectors())
#			if conf['ensemble'] == 'NPT':
#				ref_bar = float(conf['ref_bar'])
#				r_bar_freq = int(conf['ref_bar_freq'])
#				system.addForce(mm.MonteCarloBarostat(ref_bar *unit.bar, TEMP * unit.kelvin, r_bar_freq))
#			simulation = app.Simulation(topology, system, integrator, plat, prop)
#			simulation.context.setPositions(xyz)
#			simulation.context.computeVirtualSites()
#			simulation.context.applyConstraints(1e-6) # for test 
#			simulation.context.computeVirtualSites()
#			print('VECTORS: {}'.format(simulation.context.getState(getPositions=True).getPeriodicBoxVectors()))
#		else:
#			if conf['ensemble'] == 'NPT':
#				ref_bar = float(conf['ref_bar'])
#				r_bar_freq = int(conf['ref_bar_freq'])
#				system.addForce(mm.MonteCarloBarostat(ref_bar *unit.bar, TEMP * unit.kelvin, r_bar_freq))
#			system.setDefaultPeriodicBoxVectors(*vectors)
#			simulation = app.Simulation(topology, system, integrator, plat, prop)
#			simulation.context.setPositions(xyz)
#			print('VECTORS: {}'.format(simulation.context.getState(getPositions=True).getPeriodicBoxVectors()))
#	
#	if conf['md'] == 'restart':
#
#		# set up rigid bodies
#		if conf['setrigid'] != 0:
#			tmp_simu = app.Simulation(topology, system, integrator, plat, prop)
#			with open(conf['rstfile_p'],'rb') as f:
#				tmp_simu.loadCheckpoint(f)
#
#			state = tmp_simu.context.getState(getPositions=True, getVelocities=True)
#			old_vectors = state.getPeriodicBoxVectors()
#			old_positions = state.getPositions()
#			old_velocities = state.getVelocities()
#
#			v1 = old_vectors[0].value_in_unit(unit.nanometer)
#			v2 = old_vectors[1].value_in_unit(unit.nanometer)
#			v3 = old_vectors[2].value_in_unit(unit.nanometer)
#			aligned_vectors = (
#                mm.Vec3(unit.norm(v1), 0, 0),
#                mm.Vec3(v2[0], v2[1], 0),
#                mm.Vec3(v3[0], v3[1], v3[2])
#            ) * unit.nanometers
#
#			del tmp_simu
#			print('Applying constraints!')
#			system = _addrigid(system, old_positions, sys_info)
#			_integrator = mm.LangevinMiddleIntegrator(TEMP*unit.kelvin, 0.01/unit.picosecond, ts*unit.picoseconds)
#			if conf['ensemble'] == 'NPT':
#				ref_bar = float(conf['ref_bar'])
#				r_bar_freq = int(conf['ref_bar_freq'])
#				system.addForce(mm.MonteCarloBarostat(ref_bar *unit.bar, TEMP * unit.kelvin, r_bar_freq))
#			simulation = app.Simulation(topology, system, _integrator, plat, prop)
#			simulation.context.setPeriodicBoxVectors(*aligned_vectors)
#			simulation.context.setPositions(old_positions)
#			simulation.context.computeVirtualSites()
#			simulation.context.applyConstraints(1e-6)
#			simulation.context.setVelocities(old_velocities)
#		else:
#			simulation = app.Simulation(topology, system, integrator, plat, prop)
#			with open(conf['rstfile_p'],'rb') as f:
#				simulation.loadCheckpoint(f)
#
#			if conf['ensemble'] == 'NPT':
#				ref_bar = float(conf['ref_bar'])
#				r_bar_freq = int(conf['ref_bar_freq'])
#				system.addForce(mm.MonteCarloBarostat(ref_bar *unit.bar, TEMP * unit.kelvin, r_bar_freq))
#
#			state = simulation.context.getState(getPositions=True, getVelocities=True)
#			old_vectors = state.getPeriodicBoxVectors()
#			old_positions = state.getPositions()
#			old_velocities = state.getVelocities()
#
#			v1 = old_vectors[0].value_in_unit(unit.nanometer)
#			v2 = old_vectors[1].value_in_unit(unit.nanometer)
#			v3 = old_vectors[2].value_in_unit(unit.nanometer)
#			aligned_vectors = (
#                mm.Vec3(unit.norm(v1), 0, 0),
#                mm.Vec3(v2[0], v2[1], 0),
#                mm.Vec3(v3[0], v3[1], v3[2])
#            ) * unit.nanometers
#			simulation.context.setPeriodicBoxVectors(*aligned_vectors)
#			simulation.context.setPositions(old_positions)
#			simulation.context.setVelocities(old_velocities)
#
#			
#
#	if conf['md'] == 'slab':
#
#		# set up rigid bodies
#		if conf['setrigid'] != 0:
#			tmp_simu = app.Simulation(topology, system, integrator, plat, prop)
#			with open(conf['rstfile_p'],'rb') as f:
#				tmp_simu.loadCheckpoint(f)
#
#			state = tmp_simu.context.getState(getPositions=True, getVelocities=True)
#			old_vectors = state.getPeriodicBoxVectors()
#			old_positions = state.getPositions()
#			old_velocities = state.getVelocities()
#
#			new_z = old_vectors[2] * conf['vecScale']
#			new_vectors = (old_vectors[0], old_vectors[1], new_z) * unit.nanometers
#			
#			#tmp_simu.context.setPeriodicBoxVectors(*new_vectors)
#			#tmp_simu.context.setPositions(old_positions)
#			#tmp_simu.context.setVelocities(old_velocities)
#
#			del tmp_simu			
#
#			#state = simu.context.getState(getPositions=True)
#			#positions = state.getPositions()
#			print('Applying constraints!')
#			system = _addrigid(system, old_positions, sys_info)
#			_integrator = mm.LangevinMiddleIntegrator(TEMP*unit.kelvin, 0.01/unit.picosecond, ts*unit.picoseconds)
#			if conf['ensemble'] == 'NPT':
#				ref_bar = float(conf['ref_bar'])
#				r_bar_freq = int(conf['ref_bar_freq'])
#				system.addForce(mm.MonteCarloBarostat(ref_bar *unit.bar, TEMP * unit.kelvin, r_bar_freq))
#
#			simulation = app.Simulation(topology, system, _integrator, plat, prop)
#			simulation.context.setPeriodicBoxVectors(*new_vectors)
#			simulation.context.setPositions(old_positions)
#			simulation.context.setVelocities(old_velocities)
#
#		else:
#			with open(conf['rstfile_p'],'rb') as f:
#				simulation.loadCheckpoint(f)
#			if conf['ensemble'] == 'NPT':
#				ref_bar = float(conf['ref_bar'])
#				r_bar_freq = int(conf['ref_bar_freq'])
#				system.addForce(mm.MonteCarloBarostat(ref_bar *unit.bar, TEMP * unit.kelvin, r_bar_freq))
#
#
#			state = simulation.context.getState(getPositions=True, getVelocities=True)
#			old_vectors = state.getPeriodicBoxVectors()
#			old_positions = state.getPositions()
#			old_velocities = state.getVelocities()
#	
#			new_z = old_vectors[2] * conf['vecScale']
#			new_vectors = (old_vectors[0], old_vectors[1], new_z) * unit.nanometers
#	
#			simulation.context.setPeriodicBoxVectors(*new_vectors)
#			simulation.context.setPositions(old_positions)
#			simulation.context.setVelocities(old_velocities)
#			
#
#
#		simulation.minimizeEnergy()
#
#
#	# section for minimization
#	if conf['md'] == 'md':
#		print('Starting minimization !')
#		s_time = time.perf_counter()
#		simulation.minimizeEnergy()
#		positions = simulation.context.getState(getPositions=True).getPositions()
#		with open('minimized_structure.pdb', 'w') as f:
#			app.PDBFile.writeFile(topology, positions, f)	
#		e_time = time.perf_counter()
#		print('Elasping: {:.3f} seconds'.format(e_time - s_time))
#
#	if conf['md'] == 'md': # initializing the velocity according to TEMP
#		simulation.context.setVelocitiesToTemperature(TEMP * unit.kelvin)

	# add reporter for trajectory
	simulation.reporters.append(app.DCDReporter(conf['trajfile'], conf['interval'], enforcePeriodicBox=False)) # enforcePeriodicBox
																												# True = wrap; False = leave it there


	simulation.reporters.append(BoxSizeReporter('debug.log', conf['interval']))


#	simulation.reporters.append(app.StateDataReporter(
#    'debug.log', 
#    conf['interval'], 
#    step=True, 
#    potentialEnergy=True, 
#    temperature=True, 
#    volume=True,      
#    density=True,     
#    progress=True, 
#    remainingTime=True, 
#    speed=True, 
#    totalSteps=conf['step']
#))



	# chkpoint reporter/restart file
	simulation.reporters.append(app.CheckpointReporter(conf['rstfile'],reportInterval=conf['interval']))

	# energy reporter
	labels = {1:'bond',2:'angle',3:'stacking',4:'hb',5:'excluding',6:'vdw',7:'elec'}
	energies = {}
	for gid, name in labels.items():
		state = simulation.context.getState(getEnergy=True, groups={gid})
		energies[name] = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
	simulation.reporters.append(EnergyDecompositionReporter(conf['log'],reportInterval=conf['interval'],labels=labels))
	simulation.step(conf['step'])
	simulation.saveCheckpoint(conf['rstfile'])
	return None


