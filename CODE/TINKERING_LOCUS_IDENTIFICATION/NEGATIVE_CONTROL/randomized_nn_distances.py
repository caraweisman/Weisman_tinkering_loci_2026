import numpy as np
import sys
import math

gff = np.genfromtxt('randomized_dispersed_fragments.gff', dtype=str, delimiter='\t',comments=None) 
gff_sorted = gff[np.lexsort((gff[:, 3].astype(np.int64), gff[:, 0]))]

frag_dict = {} 

# add fragments to dictionaries storing fragments per contig (easier for later analysis)
for i in range(0, len(gff_sorted)):
    row = gff_sorted[i]
    contig = row[0]
    if contig not in frag_dict:
        frag_dict[contig] = []
    start = int(row[3])
    end = int(row[4])
    frag_dict[contig].append([start, end])

# make pairwise distances
pairwise_distances = [] 
for contig in frag_dict: 
    frags = frag_dict[contig]
    for i in range(0, len(frags)-1): 
        currfrag = frags[i]
        nextfrag = frags[i+1]
        pairwise_distance = nextfrag[0] - currfrag[1] 
        pairwise_distances.append(pairwise_distance)


# a) 
outfile = open('randomized_pairwise_distances', 'w') 
for dist in pairwise_distances: 
    outfile.write(str(dist) + '\n')
outfile.close()

