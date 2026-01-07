#!/bin/python

import numpy as np

def loaddump(filename):
	xyz = []
	with open(filename,'r') as f:
		while True:
			line = f.readline()
			if not line:
				break
			else:
				if len(line.split()) == 8:
					xyz.append([float(xyz) for xyz in line.split()[2:5]])

	return np.array(xyz)

def dumpxyz(xyzs, template_f, output_f):
	output = open(output_f,'w')
	with open(template_f, 'r') as f:
		while True:
			line = f.readline()
			if not line:
				break
			else:
				if len(line.split())==7:
					index = int(line.split()[0])
					chain = int(line.split()[1])
					atype = int(line.split()[2])
					q     = float(line.split()[3])
					xyz   = xyzs[index-1,:]
					print('{:3d} {:3d} {:3d} {:.3f} {:.3f} {:.3f} {:.3f}'.format(index, chain, atype, q,xyz[0], xyz[1],xyz[2]),file=output)
				else:
					print(line.strip('\n'),file=output)
					




if __name__ == '__main__':

	for i in range(0,10000001,100000):

		print(i)
		fname = 'dna.dumplammps.{}'.format(i)
		xyzs = loaddump(fname)
		print(fname,xyzs)

		dumpxyz(xyzs,'in.data','in.data.{}'.format(i))
