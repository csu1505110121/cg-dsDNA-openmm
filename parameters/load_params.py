#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

def readParam(filename):
	df = pd.read_csv(filename, sep=",", header=0)
	return df


if __name__ == '__main__':
	filename = 'protein/hps_urry_parameters.csv'
	params_hps = readParam(filename)
	print(params_hps)

	fn_dna = 'dsDNA/hps_dsDNA.csv'
	params_dna = readParam(fn_dna)
	print(params_dna)
	
