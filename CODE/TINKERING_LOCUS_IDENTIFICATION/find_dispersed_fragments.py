import numpy as np
import sys
import math
import random
import os 
import re

# usage: python nndisatances.py (main plus CLUSTERED frag gff) (contig lengths) (clusters) 

## this version filters on fragments being different than the flanking conserved genes
## but has no notion of proximity and does not output a list of fragments in a proximity range
## simpler code that will just 
## a) filter and produce a list of pairwise distances for use in the hyperexponential fitting algorithm to get initial estimates of
## insertion frequencies for "hot" and "cold" regions
## b) print out all of these filtered fragments as a gff to feed into the posterior decoding algorithm
## c) simulates random distribution of filtered and unfiltered fragments
## edit - considers whether frag is in CLUSTER with main AND does cleaner loop pass that doesn't get thrown off by overlaping conserved annotations  

genome_length = 183045933 - 20137550 # genome length minus coding region length; hard-coded. for simulation of pairwise distances

gff = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t',comments=None) # main and fragments only; code checks, but for speed
gff_sorted = gff[np.lexsort((gff[:, 3].astype(np.int64), gff[:, 0]))]

contig_lengths_file = np.genfromtxt(sys.argv[2], dtype=str)
contig_lengths = {}
for contig in contig_lengths_file: 
    contigname = contig[0]
    contiglen = int(contig[1])
    contig_lengths[contigname] = contiglen

cluster_dict = {}
with open(sys.argv[3]) as f:
    for line in f:
        fields = line.strip().split('\t')
        focal = fields[0]
        hits = fields[1:]
        cluster_dict[focal] = hits

def extract_acc(desc):
    accession = re.split('-', desc)[2]
    return accession 

# for simulation of pairwise distances: total free genome length and how to convert raw distance to 
# contig distance based on lengths of contigs
breakpoints = []
cumlen = 0
for contig in contig_lengths: 
    cumlen += contig_lengths[contig]
    breakpoints.append(cumlen) 

main_colors = ['#50C878', '#0BDA51', '#40E0D0', '#0096FF', '#5D3FD3', '#CF9FFF']
frag_colors = ['#FF2400','#E97451', '#FFAC1C','#FDDA0D']

def extract_name(desc):
    ID=re.split(';', desc)[0]
    intacc = re.split('[-_]', ID)[4:-2]
    if 'isoform' in intacc: 
        isoidx = intacc.index('isoform')
        intacc.pop(isoidx+1)
        intacc.pop(isoidx)
        intacc[isoidx-1] = intacc[isoidx-1][:-1]
    name = '-'.join(intacc)
    if 'None' in name: 
        print(desc)
    return name

contig_col = 0
start_col = 3 # start position on contig
end_col = 4 # end position on contig

# first, make a dictionary of the positions of all conserved exons 
conserved_dict = {}
for row in gff_sorted: 
    if any(s in row[8] for s in main_colors):
        name = extract_name(row[8])
        acc = extract_acc(row[8])
        contig = row[0]
        start = int(row[3])
        end = int(row[4])
        if contig in conserved_dict: 
            conserved_dict[contig].append([name, start, end, acc])
        else: 
            conserved_dict[contig] = [[name, start, end, acc]]
            
# then look through all the fragments and record their pairwise distances
# record all of their lengths for simulation later 
# make a dictionary of fragments excluding those that are plausibly derived from neighboring conserved exons
# will then go back and find pairwise distances between them
# initialize with empty vectors for every contig for simplicity, since below loop is messy

frag_dict = {} # filtered to exclude those derived from flanking conserved gene (confusing because not named as such)
for contig in contig_lengths: 
    frag_dict[contig] = []

unfiltered_frag_dict = {} # unfiltered 
for contig in contig_lengths: 
    unfiltered_frag_dict[contig] = []
    
# for ease in simulation, separate list of length of each fragment
filtered_frag_lengths = []
unfiltered_frag_lengths = []

# add fragments to dictionaries storing fragments per contig (easier for later analysis)
# two dictionaries: 
# "filtered": contains only fragments that are not derived from a flanking conserved gene (not as interesting) 
# "unfiltered": all
for i in range(0, len(gff_sorted)):
    row = gff_sorted[i]
    contig = row[0]
    start = int(row[3])
    end = int(row[4])
    name = extract_name(row[8])
    acc = extract_acc(row[8])
    if any(s in row[8] for s in frag_colors): # for all fragments
        # scan to find neighboring conserved fragments
        leftflank = None
        rightflank = None
        unfiltered_frag_dict[contig].append([name, start, end, row]) # append frag info for all fragments to unfiltered dict
        fraglen = end - start
        unfiltered_frag_lengths.append(fraglen)
        if contig in conserved_dict: # if any conserved exons are on this contig
            conserved = conserved_dict[contig] # list
            leftacc = None
            rightacc = None
            left_end = -1
            for ce in conserved:
                if ce[2] < start and ce[2] > left_end:
                    left_end = ce[2]
                    leftflank = ce[0]
                    leftacc = ce[3]
            right_start = -1
            for ce in conserved:
                if ce[1] > end:
                    if right_start == -1 or ce[1] < right_start:
                        right_start = ce[1]
                        rightflank = ce[0]
                        rightacc = ce[3]
            keep = True
            if leftflank is not None:
                if name == leftflank:
                    keep = False
                if leftacc in cluster_dict and acc in cluster_dict[leftacc]:
                    keep = False
            if rightflank is not None:
                if name == rightflank:
                    keep = False
                if rightacc in cluster_dict and acc in cluster_dict[rightacc]:
                    keep = False
            if keep:
                frag_dict[contig].append([name, start, end, row])
                filtered_frag_lengths.append(fraglen)
        else: # if no conserved exons on contig, just accept the fragment
            frag_dict[contig].append([name, start, end, row])
            
### analyses on filtered fragments
## go through and find the pairwise distances on each contig from the filtered dict for comparison to simulation
# print three output files: 
# a) pairwise distances for filtered fragments 
# b) simulated pairwise distances from random distribution of matched number/length for filtered fragments (using lengths from loop above)
# c) gff of filtered fragments (for input into superexponential fit/hmm)
filtered_pairwise_distances = [] # a
for contig in frag_dict: 
    frags = frag_dict[contig]
    ingroup = False
    for i in range(0, len(frags)-1): 
        currfrag = frags[i]
        nextfrag = frags[i+1]
        pairwise_distance = nextfrag[1] - currfrag[2] # don't append just for the filtered set 
        filtered_pairwise_distances.append(pairwise_distance)

# now, simulate distributing the filtered fragments, in their total number and length, randomly over the genome
sim_filtered_frag_raw_positions = [] # array of arrays [start, end]; will break up into contigs later
for length in filtered_frag_lengths: 
    startpos = np.random.randint(1, genome_length+1)
    endpos = startpos + length
    sim_filtered_frag_raw_positions.append([startpos, endpos])

sim_filtered_frag_raw_positions.sort(key=lambda x: x[0]) # sort in place

sim_filtered_pairwise_distances = []
for i in range(0, len(sim_filtered_frag_raw_positions)-1): 
    currfrag = sim_filtered_frag_raw_positions[i] 
    nextfrag = sim_filtered_frag_raw_positions[i+1] 
    currend = currfrag[1]
    nextstart = nextfrag[0]
    between_contigs = False
    for breakpoint in breakpoints: 
        if currend < breakpoint < nextstart: 
            between_contigs = True
    if between_contigs == False:
        pairwise_dist = nextstart - currend
        if pairwise_dist < 0: 
            sim_filtered_pairwise_distances.append(0)
        else:
            sim_filtered_pairwise_distances.append(pairwise_dist) 

## write the above a-c output files
# c)
outfile = open('dispersed_fragments.gff', 'w')
for contig in frag_dict: 
    for frag in frag_dict[contig]:
        outfile.write('\t'.join(frag[3]) + '\n')
outfile.close() 

# a) 
outfile = open('real_filtered_pairwise_distances', 'w') 
for dist in filtered_pairwise_distances: 
    outfile.write(str(dist) + '\n')
outfile.close()

# b) 
outfile = open('simulated_filtered_pairwise_distances', 'w') 
for dist in sim_filtered_pairwise_distances: 
    outfile.write(str(dist) + '\n')
outfile.close()

### analyses on UNfiltered fragments
## less than above (for now) - more interested in filtered fragments
## go through and find the pairwise distances on each contig from the filtered dict for comparison to simulation
# print two output files: 
# a) pairwise distances for unfiltered fragments 
# b) simulated pairwise distances from random distribution of matched number/length for unfiltered fragments

# for a) 
unfiltered_pairwise_frag_distances = []
for contig in unfiltered_frag_dict: 
    frags = unfiltered_frag_dict[contig]
    ingroup = False
    for i in range(0, len(frags)-1): 
        currfrag = frags[i]
        nextfrag = frags[i+1]
        pairwise_distance = nextfrag[1] - currfrag[2]
        unfiltered_pairwise_frag_distances.append(pairwise_distance)

# for b)
sim_unfiltered_frag_raw_positions = [] # array of arrays [start, end]; will break up into contigs later
for length in unfiltered_frag_lengths: 
    startpos = np.random.randint(1, genome_length+1)
    endpos = startpos + length
    sim_unfiltered_frag_raw_positions.append([startpos, endpos])

sim_unfiltered_frag_raw_positions.sort(key=lambda x: x[0]) # sort in place

sim_unfiltered_pairwise_distances = []
for i in range(0, len(sim_unfiltered_frag_raw_positions)-1): 
    currfrag = sim_unfiltered_frag_raw_positions[i] 
    nextfrag = sim_unfiltered_frag_raw_positions[i+1] 
    currend = currfrag[1]
    nextstart = nextfrag[0]
    between_contigs = False
    for breakpoint in breakpoints: 
        if currend < breakpoint < nextstart: 
            between_contigs = True
    if between_contigs == False:
        pairwise_dist = nextstart - currend
        if pairwise_dist < 0: 
            sim_unfiltered_pairwise_distances.append(0)
        else:
            sim_unfiltered_pairwise_distances.append(pairwise_dist) 
            
# write output files

# a)
outfile = open('real_unfiltered_pairwise_distances', 'w') 
for dist in unfiltered_pairwise_frag_distances: 
    outfile.write(str(dist) + '\n')
outfile.close()

# b) 
outfile = open('simulated_unfiltered_pairwise_distances', 'w') 
for dist in sim_unfiltered_pairwise_distances: 
    outfile.write(str(dist) + '\n')
outfile.close()

total_filt_frags = sum(len(v) for v in frag_dict.values())
total_unfilt_frags = sum(len(v) for v in unfiltered_frag_dict.values())

print('total filtered fragments:', total_filt_frags)
print('total unfiltered fragments:', total_unfilt_frags)
