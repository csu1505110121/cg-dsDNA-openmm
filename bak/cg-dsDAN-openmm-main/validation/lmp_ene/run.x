#!/bin/bash

for ((i=0; i<=10000000; i+=100000))
do
	echo `printf %d $i`
	python dsDNA.py 0 log.`printf %d ${i}` dna.dumplammps.`printf %d $i`
done
