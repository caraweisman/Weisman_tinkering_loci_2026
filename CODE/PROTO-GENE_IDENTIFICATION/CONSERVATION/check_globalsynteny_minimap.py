## ROW SKIPPING BUG ENSURED AGAINST
# # edit 5/17/26 to loosen flank requirement to 3000 nt on only one side

import numpy as np
import sys
import re
import warnings
warnings.filterwarnings('ignore')

# usage: python check_synteny_blast.py (minimap output file) (esl coords file) (id list for global synteny hits) (species)
# esl coord file shows the interval extracted for the synteny search and the coordinates of the fragment itself
# comparison allows determination of the relative location of the fragment in the interval, necessary for checking whether it's in the syntenic location

mmout = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t',invalid_raise=False, filling_values='', usecols=range(17))
eslcoords = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t')
idlist = np.genfromtxt(sys.argv[3], dtype=str, delimiter='\t') # necessary here but not in first step because we are doing a subset of all seqs
species = sys.argv[4]

flankseq_thresh = 3000 # longer than most RT lengths

#preindex blast results for speed
mm_index = {}
for hit in mmout:
	fid = hit[0]
	if fid not in mm_index:
		mm_index[fid] = []
	mm_index[fid].append(hit)
print('indexed dictionary')

outfile = open('2globalsynteny_results_' + species, 'w', buffering=1)
outfile.write('#Frag_ID \t Synteny_true \n')

counter = 0
for frag in eslcoords:
	fragid = frag[0]
	if fragid in idlist:
		counter += 1
		if counter%100 == 0:
			print(counter)
		frag_found = False
		coords = re.split('_', fragid)[1]
		startstop = re.split(':', coords)[1]
		start, stop = map(int, re.split('-', startstop))
		regionstart = int(frag[1])
		regionstop = int(frag[2])
		regionlen = regionstop - regionstart
		relstart = start-regionstart + 1 # minimap is 0-indexed
		relstop = stop-regionstart + 1  # minimap is 0-indexed
		if relstart - flankseq_thresh < 0: 
			leftbound = 0
		else: 
			leftbound = relstart - flankseq_thresh
		if stop + flankseq_thresh > regionstop: 
			rightbound = regionlen
		else: 
			rightbound = relstop + flankseq_thresh
		mmhits = mm_index.get(fragid, [])
		for hit in mmhits:
			hit_relstart = int(hit[2]) # left bound of alignment
			hit_relstop = int(hit[3]) # right bound of alignment
			if hit_relstart <= relstart and relstop <= hit_relstop and (hit_relstart <= leftbound or rightbound <= hit_relstop): 
				# if fragment covered by hit and flank of contiguous nt on at least one side 
				frag_found = True
		outfile.write(fragid + '\t' + str(frag_found) + '\n')
outfile.close()
