#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys, os
from pathlib import Path

current_file = Path(__file__).resolve()
project_root = current_file.parent.parent

if str(project_root) not in sys.path:
	sys.path.append(str(project_root))




def parseConf(fn,summary=False, ofn='mdp'):
	info = {
			# structure loading
			'xyz':[],
			'protein':[],
			'DNA':[],
			'fn_protein':[],
			'fn_DNA':[],
			# md parameters
			'interval':[],
			'timestep':[],
			'temperature':[],
			'dielectric':[80],
			'cellVector1':[],
			'cellVector2':[],
			'cellVector3':[],
			'cellOrigin':[],
			# machine/platform specify
			'platform':['CUDA'],
			'prop':{'Threads':'1'},
			# output log file
			'log':'log.openmm',
			'trajfile':'trajtest.dcd',
			'rstfile': 'checkpoint.chk'
			}
	with open(fn,'r') as f:
		while True:
			line = f.readline()
			if not line:
				break
			else:
				if '#' not in line and len(line.strip()):
					_record = line.split(';')[0]
					if len(_record.split()):
						record  = _record.split('=')
						args, vals = record[0].strip(), record[1].strip()
						try:
							if '.' in vals:
								vals = float(vals)
							else:
								vals = int(vals)
						except ValueError:
							pass
						info[args] = vals
	# convert to correct platform and corresponding properties			
	if info['platform']	 == 'CUDA':
		if 'DeviceIndex' in info.keys() and len(str(info['DeviceIndex'])):
			deviceidx = str(info['DeviceIndex'])
			info.pop('DeviceIndex')
		else:
			deviceidx = '0' # a default value
		if 'Precision' in info.keys() and len(info['Precision']):
			precision = info['Precision']
			info.pop('Precision')
		else:
			precision = 'mixed'
		info['prop'] = {'DeviceIndex':deviceidx, 'Precision': precision}

	if info['platform'] == 'CPU':
		if 'Threads' in info.keys() and len(str(info['Threads'])):
			threads = str(info['Threads'])
			info.pop('Threads')
		else:
			threads = '1' # set default value to 1

		info['prop'] = {'Threads':threads}

	if len(info['cellVector1']) and len(info['cellVector2']) and len(info['cellVector3']):
		v1 = [float(x) for x in info['cellVector1'].split(',')]
		v2 = [float(x) for x in info['cellVector2'].split(',')]
		v3 = [float(x) for x in info['cellVector3'].split(',')]
		info['pbc'] = ((v1[0],v1[1],v1[2]),(v2[0],v2[1],v2[2]),(v3[0],v3[1],v3[2]))
		info.pop('cellVector1')
		info.pop('cellVector2')
		info.pop('cellVector3')

	if summary:
		with open(ofn, 'w') as f:
			for k,v in info.items():
				f.write('{:10s} \t = {}\n'.format(k,v))
		

	return info


if __name__ == '__main__':
	fn_template = os.path.join(project_root,'template/md.conf')
	info = parseConf(fn_template)
	print(info)
