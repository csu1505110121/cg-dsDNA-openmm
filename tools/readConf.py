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




def parseConf(fn):
	info = {'xyz':[],
			'interval':[],
			'timestep':[],
			'temperature':[],
			'dielectric':[80],
			'cellVector1':[],
			'cellVector2':[],
			'cellVector3':[],
			'cellOrigin':[],
			'platform':['CUDA'],
			'DeviceIndex':['0'],
			'log':'log.openmm'
			}
	with open(fn,'r') as f:
		while True:
			if 
	

	return info


if __name__ == '__main__':
	fn_template = 'template/md.conf'
	info = parseConf(fn_template)
