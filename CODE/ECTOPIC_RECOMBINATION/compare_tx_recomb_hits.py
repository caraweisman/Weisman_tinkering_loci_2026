import numpy as np 
import sys

#TX_hits = np.genfromtxt('TX_hits', dtype=str, delimiter='\t', comments=None)
TX_hits = np.genfromtxt('full_TX_hits', dtype=str, delimiter='\t', comments=None)
recomb_hits = np.genfromtxt('recombination_hits', dtype=str, delimiter='\t', comments=None, skip_header=True)
no_tx_found_info = np.genfromtxt('no_tx_found', dtype=str, delimiter='\t', comments=None)
TLs = np.genfromtxt('TL_frags.gff', dtype=str, delimiter='\t', comments=None)

TL_info = {} # key = frag, value = array of results: recomb_hit, TX_hit, no tx found status 
for row in TLs:
	frag = row[8]
	TL_info[frag] = ['False', 'False', 'False'] # iniitalize all to false

for row in TX_hits: 
	frag = row[0]
	TL_info[frag][1] = 'True'

for row in recomb_hits: 
        frag = row[0]
        TL_info[frag][0] = 'True'

for frag in no_tx_found_info: 
	TL_info[frag][2] = 'True' 

# info vector components: 
# 1. is there a hit to the parental exon's flanking region?
# 2. is there a hit to the transcript?  
# 3. was there even a transcript found? makes answer to 2 irrelevant if no. (TRUE if transcript is MISSING - confusing) 

outfile = open('FULLTX_Compiled_recomb_TX_info', 'w')
#outfile = open('Compiled_recomb_TX_info', 'w')
outfile.write('#Fragment \t Recomb_hit? \t TX_hit? \t Missing_TX? \t STATUS \n')
for TL in TL_info: 
	info = TL_info[TL]
	outfile.write(TL + '\t' + info[0] + '\t' + info[1] + '\t' + info[2] + '\t') 
	if info[0] == 'True' and info[1] == 'False' and info[2] ==  'False':
	# recombination wins if there is a hit to the parental exon flanking region, but no hit to transcript, and a transcript was found
		status = 'RECOMB'
	elif info[0] == 'True' and info[1] == 'True':
		status = 'TIE'
	# tie (uninformative) if there's a hit to both parental exon flanking region and transcript
	elif info[0] == 'True' and info[1] == 'False' and info[2] ==  'True': 
		status = 'TIE'
	# tie if there's a hit to the parental exon flanking region and no transcript available for comparison
	elif info[0] == 'False' and info[1] == 'True': 
		status = 'TX'
	# TX wins if no hit to parental exon and a hit to the transcript
	elif info[0] == 'False' and info[1] == 'False' and info[2] ==  'False': 
		status = 'NONE'
	# NONE if no hit to anything
	else: 
		status = 'OTHER'
	outfile.write(status + '\n')
outfile.close()
