#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import openmm.unit as unit
from Bio.PDB.MMCIFParser import MMCIFParser


def dumpPDB(dataframe, filename):
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
		for idx, row in dataframe.iterrows():
			if len(row['bondlist']):
				f.write('{:6s}'.format('CONECT'))
				f.write('{:5d}'.format(idx+1))
				for x in row['bondlist']:
					f.write('{:5d}'.format(x+1))
				f.write('\n')


def AA2CG(fn, Mtype='protein', visual=False):
	"""
		convert all atom structure into coarse-grained ones
		
	"""
	parser = MMCIFParser()
	structure = parser.get_structure('test', fn)
	atoms = [a for a in structure.get_atoms() if a.get_name() == 'CA']
	# get the connectivity
	# 
	info = {'chain':[], 
				'name':[],   # element name
				'ATYPE':[],  # residue name
				'xyz':[],    # coordinates
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
				dumpPDB(data,'tmp.pdb')

		return data
	elif Mtype == 'DNA':
		pass
	elif Mtype == 'organic':
		pass



if __name__ == '__main__':
	data = AA2CG('../test/tdp_43.cif')
	#print(data)
