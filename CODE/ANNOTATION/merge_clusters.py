import numpy as np
import sys
from collections import defaultdict


# usage: python merge_clusters.py (extra hits file) (raw blastout file) (cluster file) 

extra_hits = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t').tolist()


distance_cutoff = 20000
overlap_thresh = 0

# make dict of blastout file for speed
col_index = 0
blastout = defaultdict(list)
with open(sys.argv[2]) as f:
	for line in f:
		hitline = line.rstrip("\n").split("\t")
		blastout[hitline[col_index]].append(hitline)


cluster_dict = {}

with open(sys.argv[3]) as f:
	for line in f:
		row = line.rstrip("\n").split("\t")
		focal_gene = row[0]
		hits = row[1:]
		cluster_dict[focal_gene] = row


count = 0
for i in range(0, len(extra_hits)-1): 
	protid = extra_hits[i][0]
	next_protid = extra_hits[i+1][0]
	startpos = min(int(extra_hits[i][6]), int(extra_hits[i][7]))
	endpos = max(int(extra_hits[i][6]), int(extra_hits[i][7]))
	next_startpos = min(int(extra_hits[i+1][6]), int(extra_hits[i+1][7]))
	next_endpos = max(int(extra_hits[i+1][6]), int(extra_hits[i+1][7]))
	contig = extra_hits[i][4]
	if protid in cluster_dict:
		distance = max(next_startpos, startpos) - min(endpos, next_endpos)
		if next_protid in cluster_dict[protid] and protid != next_protid and distance <= distance_cutoff:
			print(extra_hits[i])
			print(extra_hits[i+1])
			count = count+1
			length = int(extra_hits[i][10])
			next_length = int(extra_hits[i+1][10])
			if length >= next_length:  # use the first hit
				candhits = blastout[protid]
				for hit in candhits: 
					candstart = min(int(hit[6]), int(hit[7]))
					candend = max(int(hit[6]), int(hit[7]))
					overlap = min(next_endpos, candend) - max(next_startpos, candstart)
					if hit[4] == contig and overlap > overlap_thresh:
						extra_hits[i+1] = hit
						print('replaced')
						break
			else: 
				candhits = blastout[next_protid]
				for hit in candhits:
					candstart = min(int(hit[6]), int(hit[7]))
					candend = max(int(hit[6]), int(hit[7]))
					overlap = min(endpos, candend) - max(startpos, candstart)
					if hit[4] == contig and overlap > overlap_thresh: 
						extra_hits[i] = hit
						print('replaced')
						break



outfile = open(sys.argv[1] + '_pars', 'w')
for row in extra_hits: 
	outfile.write("\t".join(row) + "\n")
outfile.close()

print(count)
