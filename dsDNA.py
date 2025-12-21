#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from collections import Counter

import openmm as mm
import openmm.app as app
import openmm.unit as unit

from sys import stdout

class EnergyDecompositionReporter(object):
	def __init__(self, file, reportInterval, context, labels):
		self._out = open(file, 'w')
		self._interval = reportInterval
		self._context = context
		self._labels = labels
		self._out.write("Step " + " ".join(labels.values()) + "\n")

	def describeNextReport(self, simulation):
		return (self._interval, True, False, False, False)

	def report(self, simulation, state):
		step = simulation.currentStep
		self._out.write(str(step))
		for g in self._labels:
			s = self._context.getState(getEnergy=True, groups={g})
			self._out.write(f" {s.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)}")
		self._out.write("\n")

A2nm = 0.1

E_P = app.Element(number=1001, name='CGP', symbol='-p', mass=178.08)
E_A = app.Element(number=1002, name='CGA', symbol='-a', mass=134.10)
E_T = app.Element(number=1003, name='CGT', symbol='-t', mass=125.10)
E_C = app.Element(number=1004, name='CGC', symbol='-c', mass=110.10)
E_G = app.Element(number=1005, name='CGG', symbol='-g', mass=150.10)

# Make an empty topology
topology = app.Topology()


sequence1 = 'GCGTCATACAGTGC'
sequence2 = 'CGCAGTATGTCACG'


# make the chain
chain1 = topology.addChain()

# create the first residue in the chain
for seq in sequence1:
	residue = topology.addResidue(name='CG{}'.format(seq),chain=chain1)
	atom1_1 = topology.addAtom(name='CGP', element=E_P, residue=residue)
	if seq == 'A':
		atom1_2 = topology.addAtom(name='CG{}'.format(seq),element=E_A, residue=residue)
	elif seq == 'G':
		atom1_2 = topology.addAtom(name='CG{}'.format(seq),element=E_G, residue=residue)
	elif seq == 'C':
		atom1_2 = topology.addAtom(name='CG{}'.format(seq),element=E_C, residue=residue)
	elif seq == 'T':
		atom1_2 = topology.addAtom(name='CG{}'.format(seq),element=E_T, residue=residue)

	topology.addBond(atom1_1, atom1_2)
	

for i,atom in enumerate(chain1.atoms()):
	if i == 0:
		atom0 = atom
	if i !=0 and i%2 ==0 :
		topology.addBond(atom0,atom)

		atom0 = atom

# make the chain
chain2 = topology.addChain()

#print(chain2.atoms())

for seq in sequence2:
	residue = topology.addResidue(name='CG{}'.format(seq),chain=chain2)
	atom2_1 = topology.addAtom(name='CGP', element=E_P, residue=residue)
	if seq == 'A':
		atom2_2 = topology.addAtom(name='CG{}'.format(seq),element=E_A, residue=residue)
	elif seq == 'G':
		atom2_2 = topology.addAtom(name='CG{}'.format(seq),element=E_G, residue=residue)
	elif seq == 'C':
		atom2_2 = topology.addAtom(name='CG{}'.format(seq),element=E_C, residue=residue)
	elif seq == 'T':
		atom2_2 = topology.addAtom(name='CG{}'.format(seq),element=E_T, residue=residue)

	topology.addBond(atom2_1, atom2_2)
	

for i,atom in enumerate(chain2.atoms()):
	if i == 0:
		atom0 = atom
	if i !=0 and i%2 ==0 :
		topology.addBond(atom0,atom)

		atom0 = atom
#
#
#print(topology)

def readlmp(filename):
	positions = []
	with open(filename,'r') as f:
		while True:
			line = f.readline()
			if not line:
				break
			else:
				if len(line.split()) ==2 and 'atoms' in line:
					natoms = int(line.split()[0])
				if 'Atoms' in line:
					line = f.readline()				
					for i in range(natoms):
						line = f.readline()
						xyz = [float(coor) for coor in line.split()[-3:]]
						positions.append(xyz)
				if 'Bonds' in line:
					break
	return np.array(positions)
						
	
lmpdata = 'in.data'
positions = readlmp(lmpdata)				

positions = positions * A2nm * unit.nanometer

topology.setPeriodicBoxVectors(np.eye(3)*100*unit.nanometers)

with open('initial_config.pdb', 'w') as f:
	app.PDBFile.writeFile(topology, positions, f)
	

# create the system and add the particles to it
system = mm.System()
system.setDefaultPeriodicBoxVectors(*topology.getPeriodicBoxVectors())
for atom in topology.atoms():
	#print(atom)
	#print(atom.element, atom.element.mass)
	system.addParticle(atom.element.mass)

harmonic_bond_force = mm.HarmonicBondForce()

for bond in topology.bonds():
	if bond.atom1.name == 'CGP' and bond.atom2.name == 'CGP':
		#print(bond.atom1.index,bond.atom2.index, 0.55, 5000)
		harmonic_bond_force.addBond(bond.atom1.index, bond.atom2.index, 0.55*unit.nanometer, 2*50*unit.kilocalorie_per_mole/unit.angstrom**2)
		harmonic_bond_force.setUsesPeriodicBoundaryConditions(True)
	else:
		#print(bond.atom1.index,bond.atom2.index, 0.575, 5000)
		harmonic_bond_force.addBond(bond.atom1.index, bond.atom2.index, 0.575*unit.nanometer, 2*50*unit.kilocalorie_per_mole/unit.angstrom**2)
		harmonic_bond_force.setUsesPeriodicBoundaryConditions(True)

harmonic_bond_force.setForceGroup(0)


angle_expression = 'K * (cos(theta) - cos(theta0))^2'

custom_angle_force = mm.CustomAngleForce(angle_expression)
custom_angle_force.addPerAngleParameter("theta0")
custom_angle_force.addPerAngleParameter("K")

for i,atom in enumerate(chain1.atoms()):
	if atom.name == 'CGP':
		if i -2 >=0 and i +2 <len(list(chain1.atoms())):
			custom_angle_force.addAngle(atom.index-2, atom.index, atom.index+2, [180*unit.degree, 40*unit.kilocalorie_per_mole])
			custom_angle_force.setUsesPeriodicBoundaryConditions(True)
			#print('chain1',atom.index-2, atom.index, atom.index+2)


for i,atom in enumerate(chain2.atoms()):
	if atom.name == 'CGP':
		if i -2 >=0 and i +2 <len(list(chain2.atoms())):
			custom_angle_force.addAngle(atom.index-2, atom.index, atom.index+2, [180*unit.degree, 40*unit.kilocalorie_per_mole])
			custom_angle_force.setUsesPeriodicBoundaryConditions(True)
			#print('chain2',atom.index-2, atom.index, atom.index+2)

custom_angle_force.setForceGroup(1)

# section of stack interacting
stack_expression = """step(rstack-r)*
					(
					delta*(5*(sigma/r)^12-6*(sigma/r)^10)
					)
					 """
					###- delta*(5*(sigma/rstack)^12 - 6*(sigma/rstack)^10)
stack_force = mm.CustomBondForce(stack_expression)
#stack_force.addGlobalParameter('r_exclude',0.55*2**(1/6)*unit.nanometer)
stack_force.addGlobalParameter('rstack',0.62*unit.nanometer)

stack_force.addPerBondParameter("sigma")
stack_force.addPerBondParameter("delta")

for i,atom in enumerate(chain1.atoms()):
	#stack_force.addBond()
	if i % 2 ==1 and i +2 < len(list(chain1.atoms())):
		name1 = atom.name
		name2 = list(chain1.atoms())[i +2].name
		#print('chain 1',atom.index, atom.index+2)
		if name1 == name2 and name1 == 'CGA':
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 9.7*unit.kilocalorie_per_mole])
			#print('CGA-CGA',atom.index, atom.index+2)
		if name1 == name2 and name1 == 'CGT':
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 7.7*unit.kilocalorie_per_mole])
			#print('CGT-CGT',atom.index, atom.index+2)
		if name1 == name2 and name1 == 'CGC':
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 6.2*unit.kilocalorie_per_mole])
			#print('CGC-CGC',atom.index, atom.index+2)
		if name1 == name2 and name1 == 'CGG':
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 11.0*unit.kilocalorie_per_mole])
			#print('CGG-CGG',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGA','CGT']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 8.6*unit.kilocalorie_per_mole])
			#print('CGA-CGT',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGA','CGC']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 7.8*unit.kilocalorie_per_mole])
			#print('CGA-CGC',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGA','CGG']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 10.3*unit.kilocalorie_per_mole])
			#print('CGA-CGG',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGT','CGC']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 6.9*unit.kilocalorie_per_mole])
			#print('CGT-CGC',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGT','CGG']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 9.2*unit.kilocalorie_per_mole])
			#print('CGT-CGG',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGC','CGG']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 8.3*unit.kilocalorie_per_mole])
			#print('CGC-CGG',atom.index, atom.index+2)
			
		

for i,atom in enumerate(chain2.atoms()):
	#print(atom.index)
	if i % 2 ==1 and i +2 < len(list(chain2.atoms())):
		name1 = atom.name
		name2 = list(chain2.atoms())[i +2].name
		#print('chain 2',atom.index, atom.index+2)
		if name1 == name2 and name1 == 'CGA':
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 9.7*unit.kilocalorie_per_mole])
			#print('CGA-CGA',atom.index, atom.index+2)
		if name1 == name2 and name1 == 'CGT':
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 7.7*unit.kilocalorie_per_mole])
			#print('CGT-CGT',atom.index, atom.index+2)
		if name1 == name2 and name1 == 'CGC':
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 6.2*unit.kilocalorie_per_mole])
			#print('CGC-CGC',atom.index, atom.index+2)
		if name1 == name2 and name1 == 'CGG':
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 11.0*unit.kilocalorie_per_mole])
			#print('CGG-CGG',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGA','CGT']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 8.6*unit.kilocalorie_per_mole])
			#print('CGA-CGT',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGA','CGC']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 7.8*unit.kilocalorie_per_mole])
			#print('CGA-CGC',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGA','CGG']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 10.3*unit.kilocalorie_per_mole])
			#print('CGA-CGG',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGT','CGC']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 6.9*unit.kilocalorie_per_mole])
			#print('CGT-CGC',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGT','CGG']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 9.2*unit.kilocalorie_per_mole])
			#print('CGT-CGG',atom.index, atom.index+2)
		if name1 != name2 and Counter([name1, name2]) == Counter(['CGC','CGG']):
			stack_force.addBond(atom.index, atom.index+2, [0.36*unit.nanometer, 8.3*unit.kilocalorie_per_mole])
			#print('CGC-CGG',atom.index, atom.index+2)


stack_force.setUsesPeriodicBoundaryConditions(True)
stack_force.setForceGroup(2)


# section for excluded potential
excluded_expression = 'step(rcut-r)*(4*epsilon*((sigma/r)^12 - (sigma/r)^6) + epsilon); rcut=2^(1/6)*sigma'
excluded_force = mm.CustomBondForce(excluded_expression)

excluded_force.addGlobalParameter('epsilon', 4*unit.kilocalorie_per_mole)
excluded_force.addGlobalParameter('sigma', 0.55*unit.nanometer)

# add interaction between abs(i-j)==4 within the same strand
for i, atom in enumerate(chain1.atoms()):
	if atom.name != 'CGP' and i+4 < len(list(chain1.atoms())):
		excluded_force.addBond(atom.index, atom.index+4)
		print('chain 1 + Exclude force', atom.index, atom.index+4)

for i, atom in enumerate(chain2.atoms()):
	if atom.name != 'CGP' and i+4 < len(list(chain2.atoms())):
		excluded_force.addBond(atom.index, atom.index+4)
		print('chain 2 + Exclude force', atom.index, atom.index+4)

for atom1 in chain1.atoms():
	for atom2 in chain2.atoms():
		if atom1.name != 'CGP' and atom2.name != 'CGP':
			if Counter([atom1.name, atom2.name]) != Counter(['CGA','CGT']) and Counter([atom1.name, atom2.name]) != Counter(['CGC','CGG']):
				excluded_force.addBond(atom1.index,atom2.index)

# add interaction out of 1-2 and 1-3 interactions within the same strand
for i, atom1 in enumerate(chain1.atoms()):
	for j, atom2 in enumerate(chain1.atoms()):
		if i-j >4 and (atom1.name != 'CGP' and atom2.name != 'CGP'):
			if Counter([atom1.name, atom2.name]) != Counter(['CGA','CGT']) and Counter([atom1.name, atom2.name]) != Counter(['CGC','CGG']):
				excluded_force.addBond(atom1.index, atom2.index)
				print('Test chain 1 +excluded force', atom1.index,atom1.name, atom2.index, atom2.name)

for i, atom1 in enumerate(chain2.atoms()):
	for j, atom2 in enumerate(chain2.atoms()):
		if i-j >4 and (atom1.name != 'CGP' and atom2.name != 'CGP'):
			if Counter([atom1.name, atom2.name]) != Counter(['CGA','CGT']) and Counter([atom1.name, atom2.name]) != Counter(['CGC','CGG']):
				excluded_force.addBond(atom1.index, atom2.index)
				print('Test chain 1 +excluded force', atom1.index, atom1.name, atom2.index, atom2.name)


excluded_force.setUsesPeriodicBoundaryConditions(True)
excluded_force.setForceGroup(4)


# section for HB interaction

hb_expression = 'step(rhb-r)*(delta * (5*(sigma/r)^12 - 6*(sigma/r)^10))'

				####- delta*(5*(sigma/rhb)^12 - 6*(sigma/rhb)^10)
				#sigma=0.5*(sigma1+sigma2);
				#delta=sqrt(delta1*delta2)
hb_force = mm.CustomBondForce(hb_expression)
#hb_force = mm.CustomNonbondedForce(hb_expression)
#hb_force.addGlobalParameter('r_exclude',0.55*2**(1/6)*unit.nanometer)
hb_force.addGlobalParameter('rhb',0.95*unit.nanometer)

#hb_force.addPerParticleParameter('sigma')
#hb_force.addPerParticleParameter('delta')
hb_force.addPerBondParameter("sigma")
hb_force.addPerBondParameter("delta")

#for bond in topology.bonds():
#	hb_force.addExclusion(bond[0].index, bond[1].index)


#for i, atom1 in enumerate(chain1.atoms):
#	if i + 4 < len(list(chain1.atoms())):
#		if Counter([atom1.name, list(chain1.atoms())[i+4].name]) == Counter(['CGA','CGT']):

for atom1 in chain1.atoms():
	for atom2 in chain2.atoms():
		if Counter([atom1.name, atom2.name]) == Counter(['CGA','CGT']):
			hb_force.addBond(atom1.index, atom2.index, [0.60*unit.nanometer, 2.7*unit.kilocalorie_per_mole])
			#print('HB_potential',atom1.index, atom1.name, atom2.index, atom2.name)
		if Counter([atom1.name, atom2.name]) == Counter(['CGC','CGG']):
			hb_force.addBond(atom1.index, atom2.index, [0.55*unit.nanometer, 3.3*unit.kilocalorie_per_mole])
			#print('HB_potential',atom1.index, atom1.name, atom2.index, atom2.name)

for i, atom1 in enumerate(chain1.atoms()):
	for j, atom2 in enumerate(chain1.atoms()):
		if i-j >4 and Counter([atom1.name, atom2.name]) == Counter(['CGA','CGT']):
			hb_force.addBond(atom1.index, atom2.index, [0.60*unit.nanometer, 2.7*unit.kilocalorie_per_mole])
			#print('HB_potential',atom1.index, atom1.name, atom2.index, atom2.name)
		elif i-j >4 and Counter([atom1.name, atom2.name]) == Counter(['CGC','CGG']):
			hb_force.addBond(atom1.index, atom2.index, [0.55*unit.nanometer, 3.3*unit.kilocalorie_per_mole])
			#print('HB_potential',atom1.index, atom1.name, atom2.index, atom2.name)

for i, atom1 in enumerate(chain2.atoms()):
	for j, atom2 in enumerate(chain2.atoms()):
		if i-j >4 and Counter([atom1.name, atom2.name]) == Counter(['CGA','CGT']):
			hb_force.addBond(atom1.index, atom2.index, [0.60*unit.nanometer, 2.7*unit.kilocalorie_per_mole])
			#print('HB_potential',atom1.index, atom1.name, atom2.index, atom2.name)
		elif i-j >4 and Counter([atom1.name, atom2.name]) == Counter(['CGC','CGG']):
			hb_force.addBond(atom1.index, atom2.index, [0.55*unit.nanometer, 3.3*unit.kilocalorie_per_mole])
			#print('HB_potential',atom1.index, atom1.name, atom2.index, atom2.name)


hb_force.setUsesPeriodicBoundaryConditions(True)
hb_force.setForceGroup(3)

			
# section for electrostatic interaction
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
	if atom.name == 'CGP':
		elec_force.addParticle([-1])
	else:
		elec_force.addParticle([0])

# add exclusion for directly bonded atoms
elec_force.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()],1)


# set a cutoff of 3.5*ldby nm
#hb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
#hb_force.setCutoffDistance(3.5*unit.nanometers)

elec_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
elec_force.setCutoffDistance(3.5*unit.nanometer)
elec_force.setForceGroup(5)


system.addForce(harmonic_bond_force)
system.addForce(custom_angle_force)
system.addForce(stack_force)
system.addForce(hb_force)
system.addForce(excluded_force)
system.addForce(elec_force)


# run the simulation
integrator = mm.LangevinMiddleIntegrator(275*unit.kelvin, 0.01/unit.picosecond, 0.010*unit.picoseconds)
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(275*unit.kelvin)

# setup simulation reporters
# Write the trajectory to a file called 'traj.dcd'
simulation.reporters.append(app.DCDReporter('traj.dcd', 1, enforcePeriodicBox=False))

# Report information to the screen as the simulation runs
simulation.reporters.append(app.StateDataReporter(stdout, 1_000, step=True,
		potentialEnergy=True, temperature=True, volume=True, speed=True))

labels = {0:'bond',1:'angle',2:'stacking',3:'hb',4:'excluding',5:'elec'}

simulation.reporters.append(EnergyDecompositionReporter('energy_components.dat',1,simulation.context,labels))

# NVT equilibration
simulation.step(500_000)










