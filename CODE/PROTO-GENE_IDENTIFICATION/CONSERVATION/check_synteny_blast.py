import numpy as np
import sys
import re
import warnings

# usage: python check_synteny_blast.py (blast output file, self-filtered) (esl coords file) (species)
# esl coord file shows the interval extracted for the synteny search and the coordinates of the fragment itself
# comparison allows determination of the relative location of the fragment in the interval, necessary for checking whether it's in the syntenic location

blastout = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t')
eslcoords = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t')
species = sys.argv[3]

syntenic_region_matched_thresh = 0.25

#preindex blast results for speed
blast_index = {}
for hit in blastout:
	fid = hit[0]
	if fid not in blast_index:
		blast_index[fid] = []
	blast_index[fid].append(hit)
print('indexed dictionary')

outfile = open('synteny1_results_' + species, 'w', buffering=1)
outfile.write('#Frag_ID \t Syntenic_region_present_outgroup \t Syntenic_region_matched_outgroup \t Synteny_true \n')

counter = 0
for frag in eslcoords:
	counter += 1
	if counter%100 == 0:
		print(counter)
	syntenic_region_present = False 
	syntenic_region_matched = False
	synteny_found = False
	fragid = frag[0]
	coords = re.split('_', fragid)[1]
	startstop = re.split(':', coords)[1]
	start, stop = map(int, re.split('-', startstop))
	regionstart = int(frag[1])
	relstart = start-regionstart + 1 # minimap is 0-indexed
	relstop = stop-regionstart + 1  # minimap is 0-indexed
	fraglen = relstop-relstart
	blasthits = blast_index.get(fragid, [])
	for hit in blasthits:
		syntenic_region_present = True 
		if float(hit[13]) > syntenic_region_matched_thresh:
			syntenic_region_matched = True
		hit_relstart = int(hit[2])
		hit_relstop = int(hit[3])
		overlap = max(0, min(hit_relstop, relstop) - max(hit_relstart, relstart)) / (relstop - relstart) 
		if overlap >= 0.5: 
			synteny_found = True
	outfile.write(fragid + '\t' + str(syntenic_region_present) + '\t' + str(syntenic_region_matched) + '\t' + str(synteny_found) + '\n')
outfile.close()
