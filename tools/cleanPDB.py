#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from sys import stdout
import sys,os
from pathlib import Path
import argparse as ap # library for arg specified from cmd

current_file = Path(__file__).resolve()
project_root = current_file.parent.parent

if str(project_root) not in sys.path:
	sys.path.append(str(project_root))

from tools.readConf import parseConf

def _count_atoms(fn):
	atom_n_ = 0
	with open(fn, 'r') as f:
		while True:
			line = f.readline()
			if not line:
				break
			else:
				if 'ATOM' in line[:6] or 'HETATM' in line[:6]:
					atom_n_ +=1
	return atom_n_



def clean(i_conf, option='addTER'):

	# parse md configuration file
	info_conf = parseConf(i_conf)

	mix_ = info_conf['xyz']
	fnout = mix_.rsplit('.')[0]+'_{}.pdb'.format(option)
	#print(fnout)
	output = open(fnout, 'w')

	protein_ = info_conf['fn_protein']
	if protein_:
		N_prot   = info_conf['protein']
	else:
		N_prot   = 0
	dna_     = info_conf['fn_DNA']
	if dna_:
		N_DNA    = info_conf['DNA']
	else:
		N_DNA    = 0

	Mtypes = ['protein'] * int(N_prot) + ['DNA'] * int(N_DNA)
	
	if 'protein' in Mtypes:
		atom_prot_ = _count_atoms(protein_)
	if 'DNA' in Mtypes:
		atom_DNA_  = _count_atoms(dna_)

	index_ = []
	atom_count = 0
	for Mtype in Mtypes:
		if Mtype == 'protein':
			for x in range(atom_prot_):
				atom_count +=1
			index_.append(atom_count)
		elif Mtype == 'DNA':
			for x in range(atom_DNA_):
				atom_count +=1
			index_.append(atom_count)


	with open(mix_, 'r') as f:
		record_count_ = 0
		while True:
			line = f.readline()
			if not line:
				break
			else:
				if 'ATOM' in line[:6] or 'HETATM' in line[:6]:
					record_count_ +=1
					if record_count_ in index_:
						print(line.strip('\n'),file=output)
						print('TER',file=output)
					else:
						print(line.strip('\n'),file=output)



parser = ap.ArgumentParser(description='Clean structure generated')
parser.add_argument('--i_conf','-i', help='Input file of MD configuration', required=True)
args = parser.parse_args()

if __name__ == '__main__':
	try:
		clean(args.i_conf)
	except Exception as e:
		print(e)
