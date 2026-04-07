#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys, os
from sys import stdout
from pathlib import Path
import openmm as mm
import openmm.app as app
import openmm.unit as unit


current_file = Path(__file__).resolve()
project_root = current_file.parent

if str(project_root) not in sys.path:
	sys.path.append(str(project_root))

from parameters.load_params import readParam
from functions.angles import harmonic_cos_angle_form
from functions.bonds import harmonic_bond_form
from functions.nonbonded import *
from tools.io import AA2CG
from tools.model_build import build_system
from reporter.reporter import EnergyDecompositionReporter

HPS_PROT = readParam('parameters/protein/hps_urry_parameters.csv')
HPS_DNA  = readParam('parameters/dsDNA/hps_dsDNA.csv')

L = 100

#print(HPS_DNA)

data = AA2CG('test/JCTC_dsDNA.pdb',Mtype='DNA', visual=False, vis_fn='tmp.pdb')

# build the system 
topology, system = build_system(data, Mtype='DNA')

print('TEST', topology.getNumBonds())

vectors = ( (L,0,0), (0,L,0), (0,0,L) ) * unit.nanometers

system.setDefaultPeriodicBoxVectors(*vectors)
topology.setPeriodicBoxVectors(vectors)

positions = data['xyz'].values * 0.1 * unit.nanometer
#print(np.array(positions))

U_bond = harmonic_bond_form(data, Mtype='DNA')
U_angle = harmonic_cos_angle_form(data, Mtype='DNA')
U_stack = stack_12_10(data, Mtype='DNA')
U_HB	= HB_12_10(data, Mtype='DNA')
U_exclude = excluded_12_6(data, Mtype='DNA')  # have not implemented Protein
#U_stack, U_HB, U_exclude, U_vdw = stack_hb_exclude_DNA(data, Mtype='DNA')
U_elec  = elec(data, Mtype='DNA', dielectric=80)


# using those defined forces
system.addForce(U_bond)
system.addForce(U_angle)
system.addForce(U_stack)
system.addForce(U_HB)
system.addForce(U_exclude)
#system.addForce(U_vdw)
system.addForce(U_elec)


TEMP=300
plat = mm.Platform.getPlatformByName('CPU')
prop = {'Threads': '8'}

integrator = mm.LangevinMiddleIntegrator(TEMP*unit.kelvin, 0.01/unit.picosecond, 0.010*unit.picoseconds)
simulation = app.Simulation(topology, system, integrator, plat, prop)
simulation.context.setPositions(positions)

from openmm.app import PDBFile
with open('test.pdb', 'w') as f:
	PDBFile.writeFile(topology, positions, f)

#restraint_force = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
#restraint_force.addGlobalParameter("k", 100000.0 * unit.kilojoule_per_mole/unit.nanometer**2)
#restraint_force.addPerParticleParameter("x0")
#restraint_force.addPerParticleParameter("y0")
#restraint_force.addPerParticleParameter("z0")
#
#restraint_force.setForceGroup(10)
#
#for i, row in data.iterrows():
#	if row['name'] == 'PO4':
#		x, y, z = row['xyz'][0], row['xyz'][1], row['xyz'][2]
#		restraint_force.addParticle(i, [x, y, z])
#
#system.addForce(restraint_force)

#simulation.context.setParameter("k",1000)
simulation.minimizeEnergy()

positions = simulation.context.getState(getPositions=True).getPositions()


with open('minimized_structure.pdb', 'w') as f:
	PDBFile.writeFile(topology, positions, f)

#simulation.context.setParameter("k",0)
#
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

simulation.reporters.append(app.DCDReporter('traj.dcd', 10_000, enforcePeriodicBox=False))

simulation.reporters.append(app.StateDataReporter(stdout, 10_000, step=True,
   	potentialEnergy=True, temperature=True, volume=True, speed=True))


labels = {1:'bond',2:'angle',3:'stacking',4:'hb',5:'excluding',6:'vdw',7:'elec'}

energies = {}

for gid, name in labels.items():
	state = simulation.context.getState(
		getEnergy=True,
    	groups={gid}
	)
	energies[name] = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

simulation.reporters.append(EnergyDecompositionReporter('logene',reportInterval=10_000,labels=labels))

# NVT equilibration
simulation.step(1000_000)


#with open('logene','w') as f:
#	f.write('0\t')
#	for k,v in energies.items():
#		f.write('{:.6f}\t'.format(v))
#	f.write('\n')
