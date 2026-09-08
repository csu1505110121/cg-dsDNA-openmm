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
			atom_idx = idx+1
			if atom_idx > 99999:
				atom_idx = atom_idx % 99999
			f.write('{:6s}{:>5d} {:^4s} {:3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}\n'.format( \
				'ATOM', atom_idx, row['name'], row['ATYPE'], row['chain'], 1, '', \
				row['xyz'][0], row['xyz'][1], row['xyz'][2],\
				row['occ'],row['bfactor']))
	
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


def AA2CG(fn, Mtype='protein', visual=False, vis_fn='tmp.pdb', template=None, set_occ_beta=None):
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
	elif Mtype == 'nucleotide':
		atoms = [a for a in structure.get_atoms() if a.get_name() in (
									# Phosphate atoms group alpha
									"PA","O1A","O2A","O3A",
									# Phosphate atoms group beta
									"PB","O1B","O2B","O3B",
									# Phosphate atoms group gamma
									"PG","O1G","O2G","O3G",
									# sugar atoms
									"C1'","C2'","C3'","C4'","O4'",
									# ADE base
									"N1","C2","N3","C4","C5","C6","N6","N7","C8","N9",
									# THY base
									 "N1","C2","O2","N3","C4","O4","C5","C6",
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
				'anglelist':[], # angle formed with placing current bead in the middle
				'occ':[],
				'bfactor':[]
			}

	if template is not None:
		info['chain'] = template

	if set_occ_beta is not None:
		info['occ'] = set_occ_beta[Mtype]['occ'] # store occ
		info['bfactor'] = set_occ_beta[Mtype]['bfactor']   # store beta

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
			if set_occ_beta is None:
				info['occ'].append(atom.get_occupancy())  # store occ
				info['bfactor'].append(atom.get_bfactor())   # store beta

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
					# store occ and bfactor info
					if set_occ_beta is None:
						if any([a.get_occupancy() for a in atom_res if a.get_name() in PO4]):
							info['occ'].append(1)
						else:
							info['occ'].append(0)
						if any([a.get_bfactor() for a in atom_res if a.get_name() in PO4]):
							info['bfactor'].append(1)
						else:
							info['bfactor'].append(0)
					# section of BASE
					if resname in ['ADE','A']:
						if template is None:
							info['chain'].append(chain)
						info['ATYPE'].append('ADE')
						info['name'].append('ADE')
						atom_base = [a for a in atom_res if a.get_name() in ADE]
						is_occ    = any([a.get_occupancy() for a in atom_res if a.get_name() in ADE])
						is_bfactor= any([a.get_bfactor() for a in atom_res if a.get_name() in ADE])
					elif resname in ['THY','T']:
						if template is None:
							info['chain'].append(chain)
						info['ATYPE'].append('THY')
						info['name'].append('THY')
						atom_base = [a for a in atom_res if a.get_name() in THY]
						is_occ    = any([a.get_occupancy() for a in atom_res if a.get_name() in THY])
						is_bfactor= any([a.get_bfactor() for a in atom_res if a.get_name() in THY])
					elif resname in ['CYT','C']:
						if template is None:
							info['chain'].append(chain)
						info['ATYPE'].append('CYT')
						info['name'].append('CYT')
						atom_base = [a for a in atom_res if a.get_name() in CYT]
						is_occ    = any([a.get_occupancy() for a in atom_res if a.get_name() in CYT])
						is_bfactor= any([a.get_bfactor() for a in atom_res if a.get_name() in CYT])
					elif resname in ['GUA','G']:
						if template is None:
							info['chain'].append(chain)
						info['ATYPE'].append('GUA')
						info['name'].append('GUA')
						atom_base = [a for a in atom_res if a.get_name() in GUA]
						is_occ    = any([a.get_occupancy() for a in atom_res if a.get_name() in GUA])
						is_bfactor= any([a.get_bfactor() for a in atom_res if a.get_name() in GUA])
						#print(atom_base)
					com_base = com(atom_base)
					info['xyz'].append(com_base)
					info['bondlist'].append([])
					info['anglelist'].append([])
					if set_occ_beta is None:
						if is_occ:
							info['occ'].append(1)
						else:
							info['occ'].append(0)
						if is_bfactor:
							info['bfactor'].append(1)
						else:
							info['bfactor'].append(0)

		data = pd.DataFrame(data=info)


		for idx in data.index:
			pre2 = data.iloc[idx-2] if idx-2 >= 0 else None
			pre1 = data.iloc[idx-1] if idx-1 >= 0 else None
			nxt1 = data.iloc[idx+1] if idx+1 < len(data) else None
			nxt2 = data.iloc[idx+2] if idx+2 < len(data) else None

			row = data.loc[idx]

			if row['name'] == 'PO4':
				if pre2 is None and nxt2 is not None:
					data.at[idx,'bondlist'] = [idx+1,idx+2]

				elif pre2 is not None and nxt2 is not None:
					if row['chain'] == nxt2['chain'] and row['chain'] == pre2['chain']:
						data.at[idx,'bondlist'] = [idx+1, idx+2]
						data.at[idx,'anglelist'] = [idx-2, idx+2]
					elif row['chain'] == nxt2['chain'] and row['chain'] != pre2['chain']:
						data.at[idx,'bondlist'] = [idx+1,idx+2]
					else:
						data.at[idx,'bondlist'] = [idx+1]

				elif pre2 is not None and nxt2 is None:
					data.at[idx,'bondlist'] = [idx+1]

		if visual:
			dumpPDB(data,vis_fn)

		return data
	
	elif Mtype == 'nucleotide':
		chainID = None
		atomID  = 0
		chains = []
		residues = []

		ADE = ["N1","C2","N3","C4","C5","C6","N6","N7","C8","N9"]
		THY = ["N1","C2","O2","N3","C4","O4","C5","C6","C7"]
		CYT = ["N1","C2","O2","N3","C4","N4","C5","C6"]
		GUA = ["N1","C2","N2","N3","C4","C5","C6","O6","N7","C8","N9"]
		SUGAR = ["C1'","C2'","C3'","C4'","O4'"]
		PO4_A = ["PA"]
		PO4_B = ["PB"]
		PO4_G = ["PG"]

		for i, atom in enumerate(atoms):
			res = atom.get_parent()
			chain = res.get_parent()
			chains.append(chain.id)
			residues.append(res.get_id()[1])
		chains = sorted(set(chains))
		residues = set(residues)

		atoms_po4_a = [a for a in atoms if a.get_name() in PO4_A]
		atoms_po4_b = [a for a in atoms if a.get_name() in PO4_B]
		atoms_po4_g = [a for a in atoms if a.get_name() in PO4_G]
		atoms_sugar = [a for a in atoms if a.get_name() in SUGAR]
		atoms_base  = [a for a in atoms if a.get_name() in ADE+THY+CYT+GUA]

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
					# nucleotide with 3 phosphate groups, such as ATP, CTP, GTP, TTP
					if resname in ['ATP','CTP','GTP','TTP']:
						# PG
						if template is None:
							info['chain'].append(chain)
						atoms_po4_g = [a for a in atom_res if a.get_name() in PO4_G]
						com_po4_g = com(atoms_po4_g)
						info['name'].append('P2')
						info['xyz'].append(com_po4_g)
						info['ATYPE'].append(resname)
						info['bondlist'].append([1])
						info['anglelist'].append([])
						is_occ    = any([a.get_occupancy() for a in atoms_po4_g])
						is_bfactor = any([a.get_bfactor() for a in atoms_po4_g])
						if set_occ_beta is None:
							if is_occ:
								info['occ'].append(1)
							else:
								info['occ'].append(0)
							if is_bfactor:
								info['bfactor'].append(1)
							else:
								info['bfactor'].append(0)
						# PB
						if template is None:
							info['chain'].append(chain)
						atoms_po4_b = [a for a in atom_res if a.get_name() in PO4_B]
						com_po4_b = com(atoms_po4_b)
						info['name'].append('P1')
						info['xyz'].append(com_po4_b)
						info['ATYPE'].append(resname)
						info['bondlist'].append([2])
						info['anglelist'].append([0,2])
						is_occ    = any([a.get_occupancy() for a in atoms_po4_b])
						is_bfactor = any([a.get_bfactor() for a in atoms_po4_b])
						if set_occ_beta is None:
							if is_occ:
								info['occ'].append(1)
							else:
								info['occ'].append(0)
							if is_bfactor:
								info['bfactor'].append(1)
							else:
								info['bfactor'].append(0)
						# PA
						if template is None:
							info['chain'].append(chain)
						atoms_po4_a = [a for a in atom_res if a.get_name() in PO4_A]
						com_po4_a =  com(atoms_po4_a)
						info['name'].append('P1')
						info['xyz'].append(com_po4_a)
						info['ATYPE'].append(resname)
						info['bondlist'].append([3])
						info['anglelist'].append([1,3])
						is_occ    = any([a.get_occupancy() for a in atoms_po4_a])
						is_bfactor = any([a.get_bfactor() for a in atoms_po4_a])
						if set_occ_beta is None:
							if is_occ:
								info['occ'].append(1)
							else:
								info['occ'].append(0)
							if is_bfactor:
								info['bfactor'].append(1)
							else:
								info['bfactor'].append(0)
						# sugar ring
						if template is None:
							info['chain'].append(chain)
						atoms_sugar = [a for a in atom_res if a.get_name() in SUGAR]
						com_sugar = com(atoms_sugar)
						info['name'].append('CR')
						info['xyz'].append(com_sugar)
						info['ATYPE'].append(resname)
						info['bondlist'].append([4])
						info['anglelist'].append([2,4])
						is_occ    = any([a.get_occupancy() for a in atoms_sugar])
						is_bfactor = any([a.get_bfactor() for a in atoms_sugar])
						if set_occ_beta is None:
							if is_occ:
								info['occ'].append(1)
							else:
								info['occ'].append(0)
							if is_bfactor:
								info['bfactor'].append(1)
							else:
								info['bfactor'].append(0)

					# nucleotide with 2 phosphate groups, such as ADP, CDP, GDP, TDP
					if resname in ['ADP','CDP','GDP','TDP']:
						# PB
						if template is None:
							info['chain'].append(chain)
						atoms_po4_b = [a for a in atom_res if a.get_name() in PO4_B]
						com_po4_b = com(atoms_po4_b)
						info['name'].append('P2')
						info['xyz'].append(com_po4_b)
						info['ATYPE'].append(resname)
						info['bondlist'].append([1])
						info['anglelist'].append([])
						is_occ    = any([a.get_occupancy() for a in atoms_po4_b])
						is_bfactor = any([a.get_bfactor() for a in atoms_po4_b])
						if set_occ_beta is None:
							if is_occ:
								info['occ'].append(1)
							else:
								info['occ'].append(0)
							if is_bfactor:
								info['bfactor'].append(1)
							else:
								info['bfactor'].append(0)
						# PA
						if template is None:
							info['chain'].append(chain)
						atoms_po4_a = [a for a in atom_res if a.get_name() in PO4_A]
						com_po4_a =  com(atoms_po4_a)
						info['name'].append('P1')
						info['xyz'].append(com_po4_a)
						info['ATYPE'].append(resname)
						info['bondlist'].append([2])
						info['anglelist'].append([0,2])
						is_occ    = any([a.get_occupancy() for a in atoms_po4_a])
						is_bfactor = any([a.get_bfactor() for a in atoms_po4_a])
						if set_occ_beta is None:
							if is_occ:
								info['occ'].append(1)
							else:
								info['occ'].append(0)
							if is_bfactor:
								info['bfactor'].append(1)
							else:
								info['bfactor'].append(0)
						# sugar ring
						if template is None:
							info['chain'].append(chain)
						atoms_sugar = [a for a in atom_res if a.get_name() in SUGAR]
						com_sugar = com(atoms_sugar)
						info['name'].append('CR')
						info['xyz'].append(com_sugar)
						info['ATYPE'].append(resname)
						info['bondlist'].append([3])
						info['anglelist'].append([1,3])
						is_occ    = any([a.get_occupancy() for a in atoms_sugar])
						is_bfactor = any([a.get_bfactor() for a in atoms_sugar])
						if set_occ_beta is None:
							if is_occ:
								info['occ'].append(1)
							else:
								info['occ'].append(0)
							if is_bfactor:
								info['bfactor'].append(1)
							else:
								info['bfactor'].append(0)

					# nucleotide with 2 phosphate groups, such as ADP, CDP, GDP, TDP
					if resname in ['AMP','CMP','GMP','TMP']:
						# PA
						if template is None:
							info['chain'].append(chain)
						atoms_po4_a = [a for a in atom_res if a.get_name() in PO4_A]
						com_po4_a =  com(atoms_po4_a)
						info['name'].append('P2')
						info['xyz'].append(com_po4_a)
						info['ATYPE'].append(resname)
						info['bondlist'].append([1])
						info['anglelist'].append([])
						is_occ    = any([a.get_occupancy() for a in atoms_po4_a])
						is_bfactor = any([a.get_bfactor() for a in atoms_po4_a])
						if set_occ_beta is None:
							if is_occ:
								info['occ'].append(1)
							else:
								info['occ'].append(0)
							if is_bfactor:
								info['bfactor'].append(1)
							else:
								info['bfactor'].append(0)
						# sugar ring
						if template is None:
							info['chain'].append(chain)
						atoms_sugar = [a for a in atom_res if a.get_name() in SUGAR]
						com_sugar = com(atoms_sugar)
						info['name'].append('CR')
						info['xyz'].append(com_sugar)
						info['ATYPE'].append(resname)
						info['bondlist'].append([2])
						info['anglelist'].append([0,2])
						is_occ    = any([a.get_occupancy() for a in atoms_sugar])
						is_bfactor = any([a.get_bfactor() for a in atoms_sugar])
						if set_occ_beta is None:
							if is_occ:
								info['occ'].append(1)
							else:
								info['occ'].append(0)
							if is_bfactor:
								info['bfactor'].append(1)
							else:
								info['bfactor'].append(0)
						
					# BASE
					if template is None:
						info['chain'].append(chain)
					atoms_base  = [a for a in atom_res if a.get_name() in ADE+THY+CYT+GUA]
					com_base = com(atoms_base)
					if resname in ['ATP','ADP','AMP']:
						info['name'].append('XA')
						info['ATYPE'].append(resname)
						info['xyz'].append(com_base)
					elif resname in ['CTP','CDP','CMP']:
						info['name'].append('XC')
						info['ATYPE'].append(resname)
						info['xyz'].append(com_base)
					elif resname in ['GTP','GDP','GMP']:
						info['name'].append('XG')
						info['ATYPE'].append(resname)
						info['xyz'].append(com_base)
					elif resname in ['TTP','TDP','TMP']:
						info['name'].append('XT')
						info['ATYPE'].append(resname)
						info['xyz'].append(com_base)
					info['bondlist'].append([])
					info['anglelist'].append([])
					is_occ    = any([a.get_occupancy() for a in atoms_base])
					is_bfactor = any([a.get_bfactor() for a in atoms_base])
					if set_occ_beta is None:
						if is_occ:
							info['occ'].append(1)
						else:
							info['occ'].append(0)
						if is_bfactor:
							info['bfactor'].append(1)
						else:
							info['bfactor'].append(0)

		data = pd.DataFrame(data=info)

		if visual:
			dumpPDB(data,vis_fn)

		return data




def STRSUMMARY(data:pd.DataFrame) -> dict:
	"""
		summarise the system's information
		return a dictionary
	"""
	info = {'chains':[], 'atomtypes':[], 'atoms':[], 'bonds':[], 'angles':[],'total_angles':[],'rigidbodies':[]}
	#chains = []
	#atomtypes = []
	
	_rigidlist = []
	#print(data)

	for i in data.index:
		row = data.loc[i]
		Mtype = row['Mtype']
		is_rigid = row['bfactor']
		info['atoms'].append(i)
		chain_tmp = row['chain']
		if chain_tmp not in info['chains']:
			info['chains'].append(chain_tmp)
			if len(info['chains']) > 1:
				info['rigidbodies'].append(_rigidlist)
			_rigidlist = []
		elif i == len(data) -1:
			info['rigidbodies'].append(_rigidlist)
		else:
			if is_rigid:
				_rigidlist.append(i)
				
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
		elif Mtype=='nucleotide':
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

def _change_chain_name(info_conf):
	"""
		in the info_conf, the location of file are required to be specified
	"""
	chain_summary = {'protein':[],'DNA':[],'nucleotide':[]}
	if sum(info_conf['protein']):
		fn_proteins = info_conf['fn_protein']
		#print('_change_chain_name',fn_proteins)
		for fn_protein in fn_proteins:
			_data_protein = AA2CG(fn_protein,Mtype='protein')
			#print(_data_protein['chain'].values)
			#print(_data_protein)
			#_chain_protein = set(_data_protein['chain'].values)
			chain_summary['protein'].append(_data_protein['chain'].values)
		#print(chain_summary['protein'])
	if sum(info_conf['DNA']):
		#print('DNA')
		fn_DNAs	 = info_conf['fn_DNA']
		for fn_DNA in fn_DNAs:
			_data_DNA	 = AA2CG(fn_DNA,Mtype='DNA')
			#print(_data_DNA['chain'].values)
			#print(_data_DNA)
			#_chain_DNA	 = set(_data_DNA['chain'])
			chain_summary['DNA'].append(_data_DNA['chain'].values)
	if sum(info_conf['nucleotide']):
		fn_nucleotides = info_conf['fn_nucleotide']
		for fn_nucleotide in fn_nucleotides:
			_data_nucleotide = AA2CG(fn_nucleotide, Mtype='nucleotide')
			chain_summary['nucleotide'].append(_data_nucleotide['chain'].values)

	return chain_summary

def _set_occ_bfactor(info_conf):
	"""
		return the occ and bfactor from the template/monomer file
	"""
	summary = {'protein': {'occ':[],'bfactor':[]},\
				 'DNA':{'occ':[], 'bfactor':[]},\
				 'nucleotide':{'occ':[], 'bfactor':[]}}
	if sum(info_conf['protein']):
		fn_proteins = info_conf['fn_protein']
		for fn_protein in fn_proteins:
			_data_protein = AA2CG(fn_protein, Mtype='protein')
			summary['protein']['occ'].append(_data_protein['occ'].values)
			summary['protein']['bfactor'].append(_data_protein['bfactor'].values)
	if sum(info_conf['DNA']):
		fn_DNAs     = info_conf['fn_DNA']
		for fn_DNA in fn_DNAs:
			_data_DNA  = AA2CG(fn_DNA, Mtype='DNA')
			summary['DNA']['occ'].append(_data_DNA['occ'].values)
			summary['DNA']['bfactor'].append(_data_DNA['bfactor'].values)
	if sum(info_conf['nucleotide']):
		fn_nucleotides = info_conf['fn_nucleotide']
		for fn_nucleotide in fn_nucleotides:
			_data_nucleotide = AA2CG(fn_nucleotide, Mtype='nucleotide')
			summary['nucleotide']['occ'].append(_data_nucleotide['occ'].values)
			summary['nucleotide']['bfactor'].append(_data_nucleotide['bfactor'].values)

	return summary



def _splitPDB(info_conf, vis=False):
	"""
		split the assembled pdb into pdbs and assign different chains 
		and Mtype:
			- protein
			- DNA
	"""

	fn_assemble = info_conf['xyz']

	total_data = {'chain':[], 'name':[], 'ATYPE':[], 'xyz':[], 'bondlist':[], 'anglelist':[]}
	total_data = pd.DataFrame(data=total_data)

	#Mtypes = ['protein']*int(info_conf['protein']) + ['DNA']*int(info_conf['DNA']) + ['nucleotide']*int(info_conf['nucleotide'])
	Mtypes = ['protein']*sum(info_conf['protein']) + ['DNA']*sum(info_conf['DNA']) + ['nucleotide']*sum(info_conf['nucleotide'])
	Mtemplates = [i for arr in [info_conf['protein'], info_conf['DNA'], info_conf['nucleotide']] for i,n in enumerate(arr) for _ in range(n)]

	#print(Mtemplates)
	#print(info_conf)

	chain_summary = _change_chain_name(info_conf)
	_chain_protein_template = chain_summary['protein']
	_chain_DNA_template 	= chain_summary['DNA']
	_chain_nucleotide_template = chain_summary['nucleotide']

	#print('test',chain_summary)
	#print('test-protein',_chain_protein_template)

	_occ_bfactor_summary       = _set_occ_bfactor(info_conf)
	#print('_occ_bfactor_summary:',_occ_bfactor_summary)

	# estimate the total chain number
	#chains_num = info_conf['protein'] * len(set(_chain_protein_template)) + \
	#						info_conf['DNA'] * len(set(_chain_DNA_template)) + \
	#						info_conf['nucleotide'] * len(set(_chain_nucleotide_template))

	chains_num = sum([n*c for n,c in zip(info_conf['protein'], [len(set(x)) for x in _chain_protein_template])] +\
				 [n*c for n,c in zip(info_conf['DNA'], [len(set(x)) for x in _chain_DNA_template])] +\
				 [n*c for n,c in zip(info_conf['nucleotide'], [len(set(x)) for x in _chain_nucleotide_template])])
	print('Total chain number is {}'.format(chains_num))
	#print('Total chain number is {}'.format(chains_num))
	# annotate the chain number starting from 1
	chains_id = [str(i) for i in range(1, chains_num+1)]

	#print(chains_id)

	## SPLIT ASSEMBLED PDB INTO INVIDIDUAL PDBs
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
			# TER is a tag to identify the end of a certain molecule, such as protein or DNA 
			elif 'TER' in line:
				with open('_tmp.pdb','w') as tmp:
					_atm_idx = 0
					_resid = 0
					_resid_p = ''
					for record in records: # for those cases converted from gro file
						_head = record[:6]
						#_atm_idx = int(record[6:11])
						_atm_idx += 1
						_name    = record[12:16]
						_resn    = record[17:20]
						_chain   = record[21]
						_residx  = record[22:26]
						if _residx != _resid_p or _name.strip() == 'N':
							_resid_p = _residx
							_resid  += 1
						_x       = float(record[30:38])
						_y       = float(record[38:46])
						_z       = float(record[46:54])
						tmp.write('{:6s}{:>5d} {:^4s} {:3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}\n'.format(\
								_head, _atm_idx, _name, _resn, _chain, _resid, '',\
								_x, _y, _z, 1.0, 0.0))
				records = []
				Mtype = Mtypes[mol_id] # select accroding molecule type
				Mtemplate_id = Mtemplates[mol_id] # select the corresponding template id
				#print('Current template id: {}'.format(Mtemplate_id))
				#print('Current template', _chain_protein_template[Mtemplate_id] if Mtype == 'protein' else _chain_DNA_template[Mtemplate_id] if Mtype == 'DNA' else _chain_nucleotide_template[Mtemplate_id])
				_occ_bfactor = {Mtype: {'occ': _occ_bfactor_summary[Mtype]['occ'][Mtemplate_id],'bfactor': _occ_bfactor_summary[Mtype]['bfactor'][Mtemplate_id]}}
				if Mtype == 'protein':
					#print(_occ_bfactor)
					_data = AA2CG('_tmp.pdb',Mtype, visual=False, template=_chain_protein_template[Mtemplate_id], set_occ_beta=_occ_bfactor)
					#print("AA2CG is DONE")
				elif Mtype == 'DNA':
					_data = AA2CG('_tmp.pdb',Mtype, visual=False, template=_chain_DNA_template[Mtemplate_id], set_occ_beta=_occ_bfactor) 
				elif Mtype == 'nucleotide':
					_data = AA2CG('_tmp.pdb',Mtype, visual=False, template=_chain_nucleotide_template[Mtemplate_id], set_occ_beta=_occ_bfactor)

				#print('_split_data:', _data)

				_data['bondlist'] = _data['bondlist'].apply(_update_dataframe, value=len(total_data))
				_data['anglelist'] = _data['anglelist'].apply(_update_dataframe, value=len(total_data))
				_data['Mtype']	= Mtype # add another identifier `Mtype`


				if Mtype == 'protein':
					_data['chain'] = _chain_protein_template[Mtemplate_id]
				elif Mtype == 'DNA':
					_data['chain'] = _chain_DNA_template[Mtemplate_id]
				elif Mtype == 'nucleotide':
					_data['chain'] = _chain_nucleotide_template[Mtemplate_id]

				#print('Current _data\n',_data)

				unique_chains = sorted(_data['chain'].unique())
				# construct the chain map
				chain_map = {}
				for unique_chain in unique_chains:
					_x_chain = chains_id[chain_id]
					chain_id +=1
					chain_map[unique_chain]=_x_chain

				_data['chain'] = _data['chain'].map(chain_map)
					#_data['chain'] = chains_id[chain_id]

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
	#print('_split_data:',total_data)
	#total_data.to_csv('test.csv', index=False)
	return total_data
			

if __name__ == '__main__':
	# protein test
###	data = AA2CG('../test/tdp_43.cif',visual=False)
###	print(data)
###	# DNA test
#	data = AA2CG('../test/JCTC_dsDNA.pdb',Mtype='DNA', visual=True, vis_fn='tmp.pdb')
#	print(data)
###	info = STRSUMMARY(data)
###	print('chains No. {}, atomtypes No. {}'.format(len(info['chains']), len(info['atomtypes'])))
###	#for i,row in data.iterrows():
###	#	print(i, row['bondlist'] ,row['anglelist'], len(row['anglelist']))


	info_conf={'protein':0,
				'DNA': 0,
				'nucleotide': 100,
				'fn_protein':'../test/abeta.pdb',
				'fn_DNA':'../test/mix_prot_DNA/reorder.pdb',
				'fn_nucleotide':'../test/nucleotides/atp.pdb',
				'xyz':'../test/mix_ATP/mixture.pdb',
				}
#	fn_assemble = info_conf['xyz']
#	change_summary =  _change_chain_name(info_conf)
#	print(change_summary)
#	_occ_bfactor_summary = _set_occ_bfactor(info_conf)
#	print(_occ_bfactor_summary)
#	total_data = _splitPDB(info_conf,vis=True)
#	print(total_data)
#	info = STRSUMMARY(total_data)
#	print(info['chains'])
#	print(len(info['rigidbodies']))


	#AA2CG('../test/nucleotides/atp.pdb', Mtype='nucleotide', visual=True, vis_fn='tmp_ATP.pdb')
	_splitPDB(info_conf,vis=True)

#	print(total_data[total_data['Mtype']=='DNA'])
#	total_dna = total_data[total_data['Mtype']=='DNA']
#	for idx,bondlist, in zip(total_dna.index.tolist(),total_dna['bondlist'].tolist()):
#		if len(bondlist):
#			for bond in bondlist:
#				print([idx,bond],end=' ')

###	print(info['chains'])
###	print(total_data)
