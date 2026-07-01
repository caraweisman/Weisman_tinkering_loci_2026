import numpy as np 
import sys
import re
import warnings
warnings.filterwarnings("ignore")

# takes each focal TL, finds unambiguous orthologous pieces in the outgroups
# finds correspondence between them
# prints sfetch files for corresponding target/focal subregions

paf = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t', filling_values='', invalid_raise=False, usecols=range(17))  # minimap output
TL_list = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t') #sfetch file
species = sys.argv[3]

overlap_thresh = 400

#gff_focal = np.genfromtxt(sys.argv[3], dtype=str, delimiter='\t')
#gff_target = np.genfromtxt(sys.argv[4], dtype=str, delimiter='\t') 


# make a dict of all hits for each TL
TL_hit_dict = {}
for TL in TL_list: 
	TL_id = TL[0] # first col of sfetch file
	for row in paf: 
		if TL_id == row[0]: 
			if TL_id not in TL_hit_dict: 
				TL_hit_dict[TL_id] = []
			TL_hit_dict[TL_id].append(row) 



# two output files: 
# sfetch file for each subregion in the focal TL
# sfetch file for each subregion in the target TL
outfile1 = open('TL_reconstruction_esl_focal', 'w')
outfile2 = open('TL_reconstruction_esl_' + species, 'w')

# find subregions within TL (query) where there is only a single subject sequence that matches
# this is the unambiguous orthologous locus
for TL in TL_hit_dict: # TL is now the TL id, since it's the dict key from above
	hits = TL_hit_dict[TL]
	focal_contig = re.split(':', TL)[0]
	unconflicted_idx = 1
	for hit in hits: 
		unconflicted = True
		for hit2 in hits: 
			if not np.array_equal(hit, hit2):
				hitlen = int(hit[3]) - int(hit[2]) 
				hit2len = int(hit2[3]) - int(hit2[2]) 
				overlap = max(0, min(int(hit[3]), int(hit2[3])) - max(int(hit[2]), int(hit2[2])))
				if overlap > overlap_thresh or overlap >= float(hitlen)/2 or overlap >= float(hit2len)/2:
					unconflicted = False
		if unconflicted == True: # find absolute coordinates for esl-sfetch
			newid = TL + '-' + str(unconflicted_idx)
			TL_abs_coords = re.split(':', TL)[1]
			TL_abs_start = int(re.split('-', TL_abs_coords)[0])
			fochitstart = int(hit[2])
			fochitend = int(hit[3]) 
			focabshitstart = fochitstart + TL_abs_start
			focabshitend = fochitend + TL_abs_start
			# esl file for sfetch from focal species
			# -1 is for minimap off by one indexing (?)
			outfile1.write(newid + '\t' + str(focabshitstart+1) + '\t' + str(focabshitend-1) + '\t' + focal_contig + '\n')
			targethitstart = hit[7]
			targethitend = hit[8]
			targetcontig = hit[5] 
			# esl file for sfetch from target species
			# -1 is for minimap off by one indexing (?)
			outfile2.write(newid + '\t' + str(int(targethitstart)+1) + '\t' + str(int(targethitend) -1) + '\t' + targetcontig + '\n')
			unconflicted_idx +=1

outfile1.close()
outfile2.close()
