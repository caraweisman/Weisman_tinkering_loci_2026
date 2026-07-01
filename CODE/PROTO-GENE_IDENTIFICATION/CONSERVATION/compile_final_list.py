import numpy as np 
import sys
import os

species = sys.argv[1]

final_synteny_list = np.genfromtxt('final_syntenyfound_list', dtype=str, delimiter='\t') 
full_region_list = np.genfromtxt('focal_fragment_synteny_eslcoords', dtype=str, delimiter='\t', usecols=(0))

outfile = open('synteny_results_table_' + sys.argv[1], 'w', buffering=1)

for row in full_region_list: 
	if row not in final_synteny_list: 
		outfile.write(row + '\t' + str(False) + '\n') 
	else: 
		outfile.write(row + '\t' + str(True) + '\n') 

outfile.close()
