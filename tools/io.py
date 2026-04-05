#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import openmm.unit as unit
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
import re


def dumpPDB(dataframe, filename, conn=True):
	"""
		dump coarse-grained structures for visualisation	
	"""
	with open(filename,'w') as f:
		# section for xyz coordinates
		for idx, row in dataframe.iterrows():
			f.write('{:6s}{:5d} {:^4s}{:4s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n'.format( \
				'HETATM', idx+1, row['name'], row['ATYPE'], row['chain'], idx+1, '', \
				row['xyz'][0], row['xyz'][1], row['xyz'][2],\
				1.0,0.0))
	
		# section for connectivity
		if conn:
			for idx, row in dataframe.iterrows():
				if len(row['bondlist']):
					f.write('{:6s}'.format('CONECT'))
					f.write('{:5d}'.format(idx+1))
					for x in row['bondlist']:
						f.write('{:5d}'.format(x+1))
					f.write('\n')

def com(atoms):
	"""
		estimate com
	"""
	MASS = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "P": 30.974,
    "S": 32.06,
    "MG": 24.305,
    "NA": 22.990,
    "K": 39.098,
    "CA": 40.078,
    "CL": 35.45,
	}

	com = 0
	ms = 0
	for atom in atoms:
		xyz = atom.get_coord()
		m   = MASS[atom.get_name()[0]] # should be replaced with re match
		com += m*xyz
		ms +=m
	return com/ms


def AA2CG(fn, Mtype='protein', visual=False, vis_fn='tmp.pdb'):
	"""
		convert all atom structure into coarse-grained one
		current only accept two file format: PDB and CIF
		
	"""
	if fn.rsplit('.')[-1] in ['cif', 'CIF']:
		parser = MMCIFParser()
	elif fn.rsplit('.')[-1] in ['pdb', 'PDB']:
		parser = PDBParser()
	structure = parser.get_structure('test', fn)

	if Mtype == 'protein':
		atoms = [a for a in structure.get_atoms() if a.get_name() == 'CA']
	elif Mtype == 'DNA':
		atoms = [a for a in structure.get_atoms() if a.get_name() in (
									# sugar atoms
									"C1'","C2'","C3'","C4'","O4'",
									# ADE base
									"N1","C2","N3","C4","C5","C6","N6","N7","C8","N9",
									# THY base
									 "N1","C2","O2","N3","C4","O4","C5","C6","C7",
									# GUA base
									 "N1","C2","N2","N3","C4","C5","C6","O6","N7","C8","N9",
									# CYT base
									 "N1","C2","O2","N3","C4","N4","C5","C6",
																	 )]

	# get the connectivity
	# 
	info = {'chain':[], 
				'name':[],   # element name
				'ATYPE':[],  # residue name
				'xyz':[],	# coordinates
				'bondlist':[], # bond formed with current bead
				'anglelist':[] # angle formed with placing current bead in the middle
			}
	if Mtype == 'protein':
		chainID = None
		atomID  = 0
		for i,atom in enumerate(atoms):
			res = atom.get_parent()
			chain = res.get_parent()
			# if chain ID changed, the connectivity was broken and should count as new structure
			info['chain'].append(chain.id )
			info['name'].append(atom.get_name() )
			info['ATYPE'].append(res.get_resname() )
			info['xyz'].append(atom.get_coord() )
			
			nxt = atoms[i+1] if i+1 < len(atoms) else None
			pre = atoms[i-1] if i-1 >= 0 else None

			# current chain ID
			chainID = chain.id
			# next chain ID
			if nxt is not None:
				nxt_res = nxt.get_parent()
				nxt_chain = nxt_res.get_parent()
				nxt_chainID = nxt_chain.id
			# prev chain ID
			if pre is not None:
				pre_res = pre.get_parent()
				pre_chain = pre_res.get_parent()
				pre_chainID = pre_chain.id
			
			# should not count redundant info
			if pre is None and nxt is not None:
				info['anglelist'].append([])
				if chainID == nxt_chainID:
					info['bondlist'].append([atomID+1])
				else:
					info['bondlist'].append([])
			elif pre is not None and nxt is not None:
				if pre_chainID != chainID and nxt_chainID == chainID:
					info['bondlist'].append([atomID+1])
					info['anglelist'].append([])
				elif pre_chainID == chainID and nxt_chainID == chainID:
					#info['bondlist'].append([atomID-1,atomID+1])
					info['bondlist'].append([atomID+1]) # remove the redundant
					info['anglelist'].append([atomID-1,atomID+1]) # did not count pre one in bond
				elif pre_chainID == chainID and nxt_chainID != chainID:
					#info['bondlist'].append([atomID-1])
					info['bondlist'].append([]) # remove the redundant
					info['anglelist'].append([])
			elif pre is not None and nxt is None:
				if pre_chainID != chainID:
					info['bondlist'].append([])
					info['anglelist'].append([])
				elif pre_chainID == chainID:
					#info['bondlist'].append([atomID-1])
					info['bondlist'].append([]) # reduce the redundant
					info['anglelist'].append([])
			else:
				info['bondlist'].append([])
				info['anglelist'].append([])

			atomID +=1

		data = pd.DataFrame(data=info)
			#print(chain.id, res.get_resname(), res.get_id()[1], atom.get_name(),atom.get_coord())

		if visual:
			dumpPDB(data,vis_fn)

		return data


	# section for handling DNA structure
	elif Mtype == 'DNA':
		chainID = None
		atomID  = 0
		chains = []
		residues = []

		# all atom name for PO4, ADE, THY, CYT, GUA
		PO4 = ["C1'","C2'","C3'","C4'","O4'"]
		ADE = ["N1","C2","N3","C4","C5","C6","N6","N7","C8","N9"]
		THY = ["N1","C2","O2","N3","C4","O4","C5","C6","C7"]
		CYT = ["N1","C2","O2","N3","C4","N4","C5","C6"]
		GUA = ["N1","C2","N2","N3","C4","C5","C6","O6","N7","C8","N9"]

		# get the total chain ID and residue number
		for i, atom in enumerate(atoms):
			res = atom.get_parent()
			chain = res.get_parent()
			chains.append(chain.id)
			residues.append(res.get_id()[1])
		chains = set(chains)
		residues = set(residues)
		#print(chains, residues)

		for chain in chains:
			for residue in residues:
				atom_res = []
				for a in atoms:
					res = a.get_parent()
					chn = res.get_parent()
					if res.get_id()[1] == residue and chn.id == chain:
						atom_res.append(a)
						resname = res.get_resname().strip()

				# section of PO4
				info['chain'].append(chain)
				info['name'].append('PO4')
				if resname in ['ADE','A']:
					info['ATYPE'].append('ADE')
				elif resname in ['THY','T']:
					info['ATYPE'].append('THY')
				elif resname in ['CYT','C']:
					info['ATYPE'].append('CYT')
				elif resname in ['GUA','G']:
					info['ATYPE'].append('GUA')
					
				atoms_po4 = [a for a in atom_res if a.get_name() in PO4]
				com_po4 = com(atoms_po4)
				info['xyz'].append(com_po4)
				info['bondlist'].append([])
				info['anglelist'].append([])
				# section of BASE
				if resname in ['ADE','A']:
					info['chain'].append(chain)
					info['ATYPE'].append('ADE')
					info['name'].append('ADE')
					atom_base = [a for a in atom_res if a.get_name() in ADE]
				elif resname in ['THY','T']:
					info['chain'].append(chain)
					info['ATYPE'].append('THY')
					info['name'].append('THY')
					atom_base = [a for a in atom_res if a.get_name() in THY]
				elif resname in ['CYT','C']:
					info['chain'].append(chain)
					info['ATYPE'].append('CYT')
					info['name'].append('CYT')
					atom_base = [a for a in atom_res if a.get_name() in CYT]
				elif resname in ['GUA','G']:
					info['chain'].append(chain)
					info['ATYPE'].append('GUA')
					info['name'].append('GUA')
					atom_base = [a for a in atom_res if a.get_name() in GUA]
					#print(atom_base)
				com_base = com(atom_base)
				info['xyz'].append(com_base)
				info['bondlist'].append([])
				info['anglelist'].append([])

		data = pd.DataFrame(data=info)

		for idx,row in data.iterrows():
			
			pre2 = data.iloc[idx-2] if idx-2 >= 0 else None
			pre1 = data.iloc[idx-1] if idx-1 >= 0 else None
			nxt1 = data.iloc[idx+1] if idx+1 < len(data) else None
			nxt2 = data.iloc[idx+2] if idx+2 < len(data) else None

			if row['name'] == 'PO4':
				if pre2 is None and nxt2 is not None:
					row['bondlist'] = [idx+1,idx+2]

				elif pre2 is not None and nxt2 is not None:
					if row['chain'] == nxt2['chain']:
						row['bondlist'] = [idx+1, idx+2]
						row['anglelist'] = [idx+1, idx+2]
					else:
						row['bondlist'] = [idx+1]

				elif pre2 is not None and nxt2 is None:
					row['bondlist'] = [idx+1]


		if visual:
			dumpPDB(data,vis_fn)

		return data
	elif Mtype == 'organic':
		pass

def STRSUMMARY(data,Mtype='DNA'):
	chains = []
	atomtypes = []
	
	for i, row in data.iterrows():
		chain_tmp = row['chain']
		if chain_tmp not in chains:
			chains.append(chain_tmp)
		if Mtype=='DNA':
			name_tmp = row['name']
			if name_tmp not in atomtypes:
				atomtypes.append(name_tmp)
		elif Mtype=='protein':
			name_tmp = row['ATYPE']
			if name_tmp not in atomtypes:
				atomtypes.append(name_tmp)
		
	return chains, atomtypes

def readConfig(filename):
	pass



if __name__ == '__main__':
	# protein test
	data = AA2CG('../test/tdp_43.cif',visual=True)
	print(data)
	# DNA test
	data = AA2CG('../test/dsDNA.pdb',Mtype='DNA', visual=True, vis_fn='tmp.pdb')
	print(data)
	chains, atomtypes = STRSUMMARY(data)
	print('chains No. {}, atomtypes No. {}'.format(len(chains), len(atomtypes)))
	#for i,row in data.iterrows():
	#	print(i, row['bondlist'] ,row['anglelist'], len(row['anglelist']))
