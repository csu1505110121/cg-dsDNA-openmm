#!/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import re
import sys, os
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

current_file = Path(__file__).resolve()
project_root = current_file.parent.parent

if str(project_root) not in sys.path:
	sys.path.append(str(project_root))

from tools.io import STRSUMMARY
from parameters.load_params import readParam

#HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
#HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))
#HPS_NUCLEOTIDE = readParam(os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv'))



def genPSF(dataframe, filename, info_conf):
	"""
	generate the fake PSF file to visualize trajectories
	"""
	info = STRSUMMARY(dataframe)
	NCHAINS = len(info['chains'])
	NATOMS  = len(info['atoms'])
	NBONDS  = len(info['bonds'])
	NANGLES = len(info['total_angles'])

	if 'fn_hps_prot' in info_conf:
		file = Path(info_conf['fn_hps_prot'])
		if file.exists():
			HPS_PROT = readParam(info_conf['fn_hps_prot'])
		else:
			HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
	else:
		HPS_PROT = readParam(os.path.join(project_root, 'parameters/protein/hps_urry_parameters.csv'))
	if 'fn_hps_DNA' in info_conf:
		file = Path(info_conf['fn_hps_DNA'])
		if file.exists():
			HPS_DNA = readParam(info_conf['fn_hps_DNA'])
		else:
			HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))
	else:
		HPS_DNA  = readParam(os.path.join(project_root, 'parameters/dsDNA/hps_dsDNA.csv'))
	if 'fn_hps_nucleotide' in info_conf:
		file = Path(info_conf['fn_hps_nucleotide'])
		if file.exists():
			HPS_NUCLEOTIDE = readParam(info_conf['fn_hps_nucleotide'])
		else:
			HPS_NUCLEOTIDE = readParam(os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv'))
	else:
		HPS_NUCLEOTIDE = readParam(os.path.join(project_root, 'parameters/nucleotides/hps_nucleotides.csv'))

	mass_dna_map = HPS_DNA.set_index('atomtype')['mass'].to_dict()
	mass_prot_map= HPS_PROT.set_index('atomtype')['mass'].to_dict()
	mass_nucleotide_map = HPS_NUCLEOTIDE.set_index('atomtype')['mass'].to_dict()
	
	chg_dna_map  = HPS_DNA.set_index('atomtype')['charge']
	chg_prot_map = HPS_PROT.set_index('atomtype')['charge']
	chg_nucleotide_map = HPS_NUCLEOTIDE.set_index('atomtype')['charge']




	_chains = dataframe['chain'].tolist()
	_names = dataframe['name'].tolist()
	_resnames = dataframe['ATYPE'].tolist()
	_bondlists = dataframe['bondlist'].tolist()
	_anglelists = dataframe['anglelist'].tolist()
	_indices = dataframe.index.tolist()
	
	with open(filename,'w') as f:
		#----start of NTITLE-------#
		f.write('PSF CMAP\n')
		f.write('\n')
		f.write(' {} !NTITLE\n'.format(2))
		f.write('REMARKS original generated fake psf file\n')
		f.write('REMARKS contains MASS and CHARGE info currently\n')
		f.write('\n')
		#-----end of NTITLE--------#

		#----start of NATOM--------#
		f.write('{:>8d} !NATOM\n'.format(NATOMS))
		for _indice, _chain, _resname, _name in zip(_indices, _chains, _resnames, _names):
			#print(_indice, _chain, _resname, _name)
			# INDEX + CHAIN + CHAIN ID + RESNAME + ATOM NAME + ATOM TYPE + CHARGE + MASS + 0
			if _resname in HPS_PROT['atomtype'].values:
				_charge = chg_prot_map.get(_resname)
				_mass   = mass_prot_map.get(_resname)
			if _name in HPS_DNA['atomtype'].values:
				_charge = chg_dna_map.get(_name)
				_mass   = mass_dna_map.get(_name)
			if _name in HPS_NUCLEOTIDE['atomtype'].values:
				_charge = chg_nucleotide_map.get(_name)
				_mass   = mass_nucleotide_map.get(_name)
			f.write('{:>8d} {:<5s}{:<5d}{:3s}  {:<5s}{:<5s} {:>10.6f}      {:>8.4f}           0\n'.format(_indice+1, 'U', int(_chain), _resname, _name, _name, _charge, _mass))

		f.write('\n')
		#----end of NATOM----------#

		#----start of NBOND--------#
		f.write('{:>8d} !NBOND: bonds\n'.format(NBONDS))
		count = 0
		for _bond in info['bonds']:
			if len(_bond):
				count +=1
				if count % 4 ==0:
					f.write(' {:10d} {:10d}\n'.format(_bond[0]+1, _bond[1]+1))
				else:
					f.write(' {:10d} {:10d}'.format(_bond[0]+1, _bond[1]+1))
		if count % 4 !=0:
			f.write('\n')
		f.write('\n')
		#----end of MBOND----------#
		
		#----start of NTHETA-------#
		f.write('{:>8d} !NTHETA: angles\n'.format(NANGLES))
		count = 0
		for _angle in info['total_angles']:
			if len(_angle):
				count +=1
				if count %3 == 0:
					f.write(' {:10d} {:10d} {:10d}\n'.format(_angle[0]+1, _angle[1]+1, _angle[2]+1))
				else:
					f.write(' {:10d} {:10d} {:10d}'.format(_angle[0]+1, _angle[1]+1, _angle[2]+1))
		if count % 3 !=0:
			f.write('\n')
		f.write('\n')
		#----end of NTHETA --------#


		f.write('{:>8d} !NPHI: dihedrals\n'.format(0))
		f.write('\n')

		f.write('{:>8d} !NIMPHI: impropers\n'.format(0))
		f.write('\n')
		
		f.write('{:>8d} !NCRTERM: cross-terms\n'.format(0))
		f.write('\n')


if __name__ == '__main__':
	from tools.readConf import parseConf
	from tools.io import _splitPDB
	from reporter.reporter import  md
	
	info_conf = {
				'xyz':'/Volumes/T9-UOW/HPS_OpenMM/test/PM-ATP/PMx50_ATPx260/mixture.pdb',
				'protein':50,
				'DNA':0,
				'nucleotide':260,
				'fn_protein':'/Volumes/T9-UOW/HPS_OpenMM/test/PM-ATP/PM.pdb',
				'fn_DNA':'/Volumes/T9-UOW/Yi_Shen_Soy_glycinin/260422_HPS_systems/constraint/32bp.pdb',
				'fn_nucleotide':'/Volumes/T9-UOW/HPS_OpenMM/test/PM-ATP/atp.pdb',
				'cellVector1':(1000,0,0),
				'cellVector2':(0,1000,0),
				'cellVector3':(0,0,1000),
				'platform':'CPU',
				'prop':{'Threads':'8'},
				'log':'log.openmm',
				'timestep':0.01,
				'step':100,
				'pbc':((100,0,0),(0,100,0),(0,0,100)),
				'dielectric':80,
				'temperature':300,
				'C_salt': 1,
				'setrigid':0,
				'ensemble':'NVT',
				'md':'md',
				'trajfile':'traj.test.dcd',
				'rstfile': 'checkpoint.chk',
				'monitorfile':'monitor.txt',
				'interval':10,			
	}
	
	data = _splitPDB(info_conf,vis=False)
	#print(data)
	genPSF(data,'test.psf')
	#md(data, info_conf)

