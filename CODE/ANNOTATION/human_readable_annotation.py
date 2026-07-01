import numpy as np 
import sys
import re

blastfile = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t') 
protinfofile = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t')

outfile = open(sys.argv[1] + '_HR', 'w')

protdict = {}

for row in protinfofile: 
	shortname = re.split(' ', row)[0]
	fullname = row
	protdict[shortname] = fullname 

for row in blastfile: 
	name = row[0]
	newname = protdict[name]
	newrow = [newname] + row[1:].tolist()
	rowtowrite = '\t'.join(newrow)
	outfile.write(rowtowrite + '\n')
outfile.close()
