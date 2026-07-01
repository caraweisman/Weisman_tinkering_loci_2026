import numpy as np
import glob

prefix = 'focal_TL_coverage_'
species_list = ['Dhel', 'Dazt', 'Dalg', 'Dath']
coverage_files = [prefix + i for i in species_list]

TL_list = np.genfromtxt('TL_sfetch', dtype=str, delimiter='\t', usecols=(0))

TL_coverage_dict = {}
for TL in TL_list: 
	TL_coverage_dict[TL] = [] # list will be each species' covearge in order of specieslist

for file in coverage_files: 
	species_TL_coverage_dict = {} 
	coverage_list = np.genfromtxt(file, dtype=str, delimiter='\t')
	for row in coverage_list: 
		species_TL_coverage_dict[row[0]] = row[1] # key = TL, value = coverage
	for TL in TL_list: 
		coverage = species_TL_coverage_dict.get(TL, 0) 
		TL_coverage_dict[TL].append(coverage)

outfile = open('compiled_coverage', 'w')
outfile.write('TL' + '\t' + "\t".join(species_list) + '\n')

for TL in TL_coverage_dict:
	outfile.write(TL) 
	coverage_vals = TL_coverage_dict[TL]
	for spec in coverage_vals:
		outfile.write('\t' + str(spec)) 
	outfile.write('\n') 
outfile.close()
