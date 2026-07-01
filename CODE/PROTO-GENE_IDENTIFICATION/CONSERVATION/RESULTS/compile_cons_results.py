import glob
import numpy as np

# for new code with unfiltered

results_file_prefix = 'synteny_results_table_'

species = ['Dath', 'Dalg', 'Dazt', 'Dhel']

results_files = []
for s in species: 
	results_files.append(results_file_prefix + s)

outfile = open('conservation_results', 'w')
outfile.write('Fragment' + '\t')

result_dict = {}
for file in results_files: 
	outfile.write(file + '\t') 
	results = np.genfromtxt(file, dtype=str, delimiter='\t')
	for result in results: 
		fragid = result[0]
		state = result[1]
		if fragid not in result_dict: 
			result_dict[fragid] = []
		result_dict[fragid].append(state)
outfile.write('\n') 

consoutfile = open('all_unconserved_fraglist', 'w')

for fragid in result_dict: 
	cons = False
	results = result_dict[fragid]
	outfile.write(fragid)
	for result in results: 
		outfile.write('\t' + result)
		if result == 'True': 
			cons = True
	outfile.write('\n') 
	if cons == False: 
		consoutfile.write(fragid + '\n')
outfile.close()
consoutfile.close()

outfile = open('age_list', 'w')

for fragid in result_dict: 
	results = result_dict[fragid]
	outfile.write(fragid + '\t')
	if results[3] == 'True': 
		outfile.write('DHEL')
	elif results[2] == 'True': 
		outfile.write('DAZT')
	elif results[1] == 'True' or results[0] == 'True':
		outfile.write('DALG/DATH')
	else: 
		outfile.write('DAFF')
	outfile.write('\n') 
outfile.close()
