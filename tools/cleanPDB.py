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
	#print(info_conf)

	mix_ = info_conf['xyz']
	fnout = mix_.rsplit('.')[0]+'_{}.pdb'.format(option)
	#print(fnout)
	output = open(fnout, 'w')

	Mtypes = ['protein'] * sum(info_conf['protein']) + ['DNA'] * sum(info_conf['DNA']) + ['nucleotide'] * sum(info_conf['nucleotide'])
	Mtemplates = [i for arr in [info_conf['protein'], info_conf['DNA'], info_conf['nucleotide']] for i,n in enumerate(arr) for _ in range(n)]

	#print(Mtypes)
	#print(Mtemplates)

	atom_prots_ 		 = []
	atom_DNAs_  		 = []
	atom_nucleotides_ = []

	if sum(info_conf['protein']):
		for fn_protein in info_conf['fn_protein']:
			atom_prots_.append(_count_atoms(fn_protein))

	if sum(info_conf['DNA']):
		for fn_DNA in info_conf['fn_DNA']:
			atom_DNAs_.append(_count_atoms(fn_DNA))

	if sum(info_conf['nucleotide']):
		for fn_nucleotide in info_conf['fn_nucleotide']:
			atom_nucleotides_.append(_count_atoms(fn_nucleotide))


	index_ = []
	atom_count = 0
	for Mtype,Mtemplate in zip(Mtypes,Mtemplates):
		if Mtype == 'protein':
			atom_prot_ = atom_prots_[Mtemplate]
			for x in range(atom_prot_):
				atom_count +=1
			index_.append(atom_count)
		elif Mtype == 'DNA':
			atom_DNA_ = atom_DNAs_[Mtemplate]
			for x in range(atom_DNA_):
				atom_count +=1
			index_.append(atom_count)
		elif Mtype == 'nucleotide':
			atom_nucleotide_ = atom_nucleotides_[Mtemplate]
			for x in range(atom_nucleotides_[Mtemplate]):
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
