#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys, os
from sys import stdout
from pathlib import Path
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import argparse as ap # library for arg specified from cmd

current_file = Path(__file__).resolve()
project_root = current_file.parent

if str(project_root) not in sys.path:
	sys.path.append(str(project_root))

from parameters.load_params import readParam
from functions.angles import harmonic_cos_angle_form
from functions.bonds import harmonic_bond_form
from functions.nonbonded import *
from tools.io import AA2CG, STRSUMMARY, _splitPDB
from tools.model_build import build_system
from tools.readConf import parseConf
from tools.converting import generateDATA
from reporter.reporter import EnergyDecompositionReporter, md

HPS_PROT = readParam(os.path.join(project_root,'parameters/protein/hps_urry_parameters.csv'))
HPS_DNA  = readParam(os.path.join(project_root,'parameters/dsDNA/hps_dsDNA.csv'))


def main_exe(i_conf,dumplammps,verbose):
	# step1: parse md configuration file
	if verbose:
		print('Step 1: Parsing the input file')
	info_conf = parseConf(i_conf,summary=verbose, ofn='mdp')
	#print(info_conf)

	# step2: parse system and coarse grained them according to rules
	if verbose:
		print('Step 2: Parsing the assembled pdb')
	data      = _splitPDB(info_conf,vis=verbose)
	if dumplammps:
		generateDATA(data, info_conf)
	
	# step3: run simulation
	if verbose:
		print('Step 3: Runing simulations')
	if info_conf['md'] != 'sp':
		md(data, info_conf)
	elif info_conf['md'] == 'sp':
		_sp_estimation(data, info_conf)



parser = ap.ArgumentParser(description='Coarse-grained simulation')
parser.add_argument('--i_conf','-i', help='Input file of MD configuration', required=True)
parser.add_argument('--dumplammps','-lammps', action='store_true', help='Dump files for LAMMPS simulation')
parser.add_argument('--verbose','-v',action='store_true', help='Print verbose information')
args = parser.parse_args()

if __name__ == '__main__':
	try:
		main_exe(args.i_conf, args.dumplammps, args.verbose)
	except Exception as e:
		print(e)
