import numpy as np
import sys
import os
import random

gff = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t') 
contiglenfile = np.genfromtxt(sys.argv[2], dtype=str)

contig_lengths = {}
for line in contiglenfile: 
	contig = line[0]
	length = int(line[1])
	contig_lengths[contig] = length


frag_lengths = []
for line in gff: 
	fragstart = int(line[3])
	fragend = int(line[4])
	fraglen = fragend - fragstart
	frag_lengths.append(fraglen)

def pick_contig():
	contigs = list(contig_lengths.keys())
	lengths = list(contig_lengths.values())
	return random.choices(contigs, weights=lengths, k=1)[0]

sim_contig_frags = {}
for contig in contig_lengths: 
	sim_contig_frags[contig] = []


outfile = open('randomized_dispersed_fragments.gff', 'w')
for fraglen in frag_lengths: 
	contig = pick_contig()
	contiglen = contig_lengths[contig]
	startpos = random.randint(1, contiglen-fraglen)
	endpos = startpos + fraglen
	outfile.write(contig + '\t' + 'tblastn' + '\t' + 'match' + '\t' + str(startpos) + '\t' + str(endpos) + '\t' + '.' + '\t' + '+' + '\t' + '.' + '\t' + 'ID=SIMULATED' + '\n') 
outfile.close()
