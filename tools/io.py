#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import openmm.unit as unit
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
import re
import sys, os
import warnings
warnings.filterwarnings("ignore")

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


def _reorder_resid(structure):
	"""
		reorder residue ID
	"""
	res_counter = 0
	last_res_obj = None

	all_atoms = list(structure.get_atoms())

	for atom in all_atoms:
		res_obj = atom.get_parent()
		if res_obj is not last_res_obj:
			res_counter +=1
			old_id = res_obj.id
			res_obj.id = (old_id[0], 10000 + res_counter, old_id[2]) # initially assign large residue to avoid conflicts
			last_res_obj = res_obj

	res_counter = 0
	last_res_obj = None
	for atom in all_atoms:
		res_obj = atom.get_parent()
		if res_obj is not last_res_obj:
			res_counter +=1
			old_id = res_obj.id
			res_obj.id = (old_id[0], res_counter, old_id[2])
			last_res_obj = res_obj

	return structure


def AA2CG(fn, Mtype='protein', visual=False, vis_fn='tmp.pdb', template=None):
	"""
		convert all atom structure into coarse-grained one
		current only accept two file format: PDB and CIF
		
	"""
	if fn.rsplit('.')[-1] in ['cif', 'CIF']:
		parser = MMCIFParser()
	elif fn.rsplit('.')[-1] in ['pdb', 'PDB']:
		parser = PDBParser()
	structure = parser.get_structure('test', fn)

	# reorder residue ID
	# this step is replaced by pdb4amber
	#structure = _reorder_resid(structure)


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

	if template is not None:
		info['chain'] = template

	if Mtype == 'protein':
		chainID = None
		atomID  = 0
		for i,atom in enumerate(atoms):
			res = atom.get_parent()
			chain = res.get_parent()
			# if chain ID changed, the connectivity was broken and should count as new structure
			if template is None:
				info['chain'].append(chain.id )
			info['name'].append(atom.get_name() )
			if res.get_resname() in ['HIS','HIE','HID']: # just convert HIE and HID to HIS
				info['ATYPE'].append('HIS')
			else:
				info['ATYPE'].append(res.get_resname())
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
		chains = sorted(set(chains))
		residues = set(residues)
		
		#print(chains)

		for chain in chains:
			for residue in residues:
				atom_res = []
				for a in atoms:
					res = a.get_parent()
					chn = res.get_parent()
					if res.get_id()[1] == residue and chn.id == chain:
						atom_res.append(a)
						resname = res.get_resname().strip()

				if len(atom_res):
					# section of PO4
					if template is None:
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
						if template is None:
							info['chain'].append(chain)
						info['ATYPE'].append('ADE')
						info['name'].append('ADE')
						atom_base = [a for a in atom_res if a.get_name() in ADE]
					elif resname in ['THY','T']:
						if template is None:
							info['chain'].append(chain)
						info['ATYPE'].append('THY')
						info['name'].append('THY')
						atom_base = [a for a in atom_res if a.get_name() in THY]
					elif resname in ['CYT','C']:
						if template is None:
							info['chain'].append(chain)
						info['ATYPE'].append('CYT')
						info['name'].append('CYT')
						atom_base = [a for a in atom_res if a.get_name() in CYT]
					elif resname in ['GUA','G']:
						if template is None:
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
					if row['chain'] == nxt2['chain'] and row['chain'] == pre2['chain']:
						row['bondlist'] = [idx+1, idx+2]
						row['anglelist'] = [idx-2, idx+2]
					elif row['chain'] == nxt2['chain'] and row['chain'] != pre2['chain']:
						row['bondlist'] = [idx+1,idx+2]
					else:
						row['bondlist'] = [idx+1]

				elif pre2 is not None and nxt2 is None:
					row['bondlist'] = [idx+1]


		if visual:
			dumpPDB(data,vis_fn)

		return data
	elif Mtype == 'organic':
		pass

def STRSUMMARY(data,Mtype='DNA') -> dict:
	"""
		summarise the system's information
		return a dictionary
	"""
	info = {'chains':[], 'atomtypes':[], 'atoms':[], 'bonds':[], 'angles':[],'total_angles':[]}
	#chains = []
	#atomtypes = []
	
	for i, row in data.iterrows():
		Mtype = row['Mtype']
		info['atoms'].append(i)
		chain_tmp = row['chain']
		if chain_tmp not in info['chains']:
			info['chains'].append(chain_tmp)
		bonds = row['bondlist']
		angles = row['anglelist']
		if len(bonds):
			for bond in bonds:
				info['bonds'].append([i,bond])
		if len(angles) and Mtype == 'DNA':  # only store angles related to 'DNA'
			info['angles'].append([angles[0], i, angles[1]])
		if len(angles): # store all angles in the system
			info['total_angles'].append([angles[0], i, angles[1]])
			
		if Mtype=='DNA':
			name_tmp = row['name']
			if name_tmp not in info['atomtypes']:
				info['atomtypes'].append(name_tmp)
		elif Mtype=='protein':
			name_tmp = row['ATYPE']
			if name_tmp not in info['atomtypes']:
				info['atomtypes'].append(name_tmp)
		
	return info


def _update_dataframe(columns, value):
	"""
		reindex the index of different molecules,
		- where columns specify the [anglelist] and [bondlist]
		- value is the length of previously accumulated atom number
	"""
	if isinstance(columns, list) and len(columns) >0:
		return [i + value for i in columns]
	return columns

def _change_chain_name(fn_assemble, info_conf):
	"""
		in the info_conf, the location of file are required to be specified
	"""
	chain_summary = {'protein':[],'DNA':[]}
	if info_conf['protein']:
		fn_protein = info_conf['fn_protein']
		_data_protein = AA2CG(fn_protein,Mtype='protein')
		#print(_data_protein['chain'].values)
		#print(_data_protein)
		#_chain_protein = set(_data_protein['chain'].values)
		chain_summary['protein'] = _data_protein['chain'].values
	if info_conf['DNA']:
		fn_DNA	 = info_conf['fn_DNA']
		_data_DNA	 = AA2CG(fn_DNA,Mtype='DNA')
		#print(_data_DNA['chain'].values)
		#print(_data_DNA)
		#_chain_DNA	 = set(_data_DNA['chain'])
		chain_summary['DNA']	 = _data_DNA['chain'].values

	return chain_summary


def _splitPDB(fn_assemble,info_conf, vis=False):
	"""
		split the assembled pdb into pdbs and assign different chains 
		and Mtype:
			- protein
			- DNA
	"""

	total_data = {'chain':[], 'name':[], 'ATYPE':[], 'xyz':[], 'bondlist':[], 'anglelist':[]}
	total_data = pd.DataFrame(data=total_data)

	Mtypes = ['protein']*int(info_conf['protein']) + ['DNA']*int(info_conf['DNA'])
	#print(Mtypes)

	chain_summary = _change_chain_name(fn_assemble, info_conf)
	_chain_protein_template = chain_summary['protein']
	_chain_DNA_template 	= chain_summary['DNA']

	#print(_chain_protein_template, _chain_DNA_template)

	# estimate the total chain number
	chains_num = info_conf['protein'] * len(set(_chain_protein_template)) + info_conf['DNA'] * len(set(_chain_DNA_template))
	#print('Total chain number is {}'.format(chains_num))
	# annotate the chain number starting from 1
	chains_id = [str(i) for i in range(1, chains_num+1)]
	
	mol_id = 0
	chain_id = 0
	with open(fn_assemble,'r') as f:
		records = []
		while True:
			line = f.readline()
			if not line:
				break
			elif 'ATOM' in line or 'HETATM' in line:
				records.append(line)
			elif 'TER' in line:
				with open('_tmp.pdb','w') as tmp:
					for record in records:
						tmp.write(record)
				records = []
				Mtype = Mtypes[mol_id] # select accroding molecule type
				if Mtype == 'protein':
					_data = AA2CG('_tmp.pdb',Mtype, visual=False, template=_chain_protein_template) 
				elif Mtype == 'DNA':
					_data = AA2CG('_tmp.pdb',Mtype, visual=False, template=_chain_DNA_template) 
					
				_data['bondlist'] = _data['bondlist'].apply(_update_dataframe, value=len(total_data))
				_data['anglelist'] = _data['anglelist'].apply(_update_dataframe, value=len(total_data))
				_data['Mtype']	= Mtype # add another identifier `Mtype`

				if Mtype == 'protein':
					_data['chain'] = _chain_protein_template
				elif Mtype == 'DNA':
					_data['chain'] = _chain_DNA_template

				unique_chains = sorted(_data['chain'].unique())
				# construct the chain map
				chain_map = {}
				for unique_chain in unique_chains:
					_x_chain = chains_id[chain_id]
					chain_id +=1
					chain_map[unique_chain]=_x_chain

				_data['chain'] = _data['chain'].map(chain_map)
					#_data['chain'] = chains_id[chain_id]
			#	elif Mtype == 'DNA':
			#		unique_chains = sorted(_data['chain'].unique())
			#		print(unique_chains)
			#		_left_chain = chains_id[chain_id]
			#		chain_id +=1
			#		_right_chain = chains_id[chain_id]
			#		chain_map = {unique_chains[0]:_left_chain, unique_chains[1]:_right_chain}
			#		_data['chain'] = _data['chain'].map(chain_map)

				mol_id +=1
				#chain_id +=1
				total_data = pd.concat([total_data, _data], ignore_index=True)
				# removing the temp file, such as _tmp.pdb
				if os.path.exists('_tmp.pdb'):
					os.remove('_tmp.pdb')
	if vis:
		total_data_new = total_data.copy()
		total_data_new['chain']	= 'A'
		dumpPDB(total_data_new, filename='total_mix.pdb')					

	return total_data
			

if __name__ == '__main__':
	# protein test
###	data = AA2CG('../test/tdp_43.cif',visual=False)
###	print(data)
###	# DNA test
###	data = AA2CG('../test/JCTC_dsDNA.pdb',Mtype='DNA', visual=True, vis_fn='tmp.pdb')
###	print(data)
###	info = STRSUMMARY(data)
###	print('chains No. {}, atomtypes No. {}'.format(len(info['chains']), len(info['atomtypes'])))
###	#for i,row in data.iterrows():
###	#	print(i, row['bondlist'] ,row['anglelist'], len(row['anglelist']))


	info_conf={'protein':100,
				'DNA': 25,
				'fn_protein':'../test/abeta.pdb',
				'fn_DNA':'../test/mix_prot_DNA/reorder.pdb'
				}
	fn_assemble = '../test/mix_prot_DNA/mixture.pdb'
	change_summary =  _change_chain_name(fn_assemble, info_conf)
	print(change_summary)
	total_data = _splitPDB(fn_assemble,info_conf,vis=True)
	info = STRSUMMARY(total_data)
	#print(info)

###	print(info['chains'])
###	print(total_data)
