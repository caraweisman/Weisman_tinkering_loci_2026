import sys
import numpy as np
from itertools import product
import re 

# usage: check_conserved.py (focal blastout - from coding analysis in above subdir to find chimeric orfs) (target blastout - from this dir) (esl-translate ORF fasta from reads to link ORF to reads) (reads to locus _withreads)

focal_blastout = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t')
esl_peps = sys.argv[2]
reads_to_locus = sys.argv[3]
target_blastout_paths = sys.argv[4:]

#target_blastout = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t')


MAX_GAP = 30000

def extract_name(desc):
    # desc looks like: "NP_001303394.1 nicotinic acetylcholine receptor alpha4, isoform H [Drosophila melanogaster]"
    # strip the leading accession (first whitespace token)
    rest = desc.split(' ', 1)[1] if ' ' in desc else desc
    # remove ", isoform X" (X is one or more non-bracket, non-comma chars)
    rest = re.sub(r',\s*isoform\s+[^,\[]+', '', rest)
    # collapse all whitespace to single dashes
    name = re.sub(r'\s+', '-', rest.strip())
    return name

orf_focal_ranges = {}
for row in focal_blastout: 
    orf = row[0]
    if orf not in orf_focal_ranges: 
        orf_focal_ranges[orf] = []
    hitstart = int(row[2]) 
    hitstop = int(row[3]) 
    hitid = row[4]
    orf_focal_ranges[orf].append([hitstart, hitstop, hitid]) 


orf_to_read_dict = {}
with open(esl_peps) as f:
    for line in f:
        fields = line.split()
        if '>' == fields[0][0]: 
            orf = fields[0][1:]
            read = re.split('=', fields[1])[1]
            orf_to_read_dict[orf] = read

reads_to_locus_dict = {} # key = read, value = locus row index
locus_dict = {} # will store supporting read counts here. key = locus row index, value = [read, orf components comma sep]
with open(reads_to_locus) as f:
    for line_idx, line in enumerate(f):
        fields = line.rstrip('\n').split('\t')
        prefix = fields[:-1]
        readlist = re.split(',', fields[-1])
        for read in readlist:
            reads_to_locus_dict[read] = line_idx
        locus_dict[line_idx] = [] # initialize

orf_to_protid_dict = {} # key = orfid, value = [prots]

# Per-orf tracking across all targets
orf_conserved_anywhere = {orf: False for orf in orf_focal_ranges}

for target_path in target_blastout_paths:
    target_blastout = np.genfromtxt(target_path, dtype=str, delimiter='\t')
    orf_target_blasthits = {}
    for row in target_blastout: 
        orf = row[0]
        orfstart = int(row[2])
        orfstop = int(row[3])
        contig = row[4]
        contigstart = min(int(row[6]), int(row[7]))
        contigstop = max(int(row[6]), int(row[7]))
        if int(row[6]) < int(row[7]): 
            strand = '+'
        else: 
            strand = '-'
        if orf not in orf_target_blasthits: 
            orf_target_blasthits[orf] = {} # format: key = orf, value = ANOTHER DICT with contig as key, [[orfstart, orfend, contigstart, contigend]] as value
        if contig not in orf_target_blasthits[orf]: 
            orf_target_blasthits[orf][contig] = []
        orf_target_blasthits[orf][contig].append([orfstart, orfstop, contigstart,contigstop, strand]) 
    for orf in orf_focal_ranges: 
        prot_ranges = sorted(orf_focal_ranges[orf], key=lambda r: r[0])  # sort prot ranges by prot start so that order can be easily assessed below
        prot_ids = [extract_name(i[2]) for i in prot_ranges]
        orf_to_protid_dict[orf] = ','.join(prot_ids)
        orf_protids = ','.join(sorted(prot_ids))
        target_blasthit_contigs = orf_target_blasthits[orf]
        successful_combos = []
        orf_prox = False
        orf_strand = False
        orf_order = False
        for contig in target_blasthit_contigs: 
            hits = target_blasthit_contigs[contig]
            prot_ranges_found = [False] * len(prot_ranges)
            candidates = [[] for _ in prot_ranges] # keep track of hits for each orf range
            for hit in hits: 
                targetorfstart, targetorfstop, contigstart, contigstop  = hit[0], hit[1], hit[2], hit[3]
                for i, (focal_start, focal_stop, _hitid) in enumerate(prot_ranges):
                    overlap = min(focal_stop, targetorfstop) - max(focal_start, targetorfstart) + 1
                    if overlap >= 40 or overlap >= 0.5 * (focal_stop - focal_start + 1):
                        prot_ranges_found[i] = True
                        candidates[i].append(hit)
            if not all(prot_ranges_found): # first just check whether every protein region is even on the contig; many will fail here
                 continue
            prox = False
            prox_combos = []
            for combo in product(*candidates): # amazing itertools thing that makes all possible products of hits per orf range so they can be tested for contiguity
                sorted_combo = sorted(combo, key=lambda h: h[2])
                if all(sorted_combo[i+1][2] - sorted_combo[i][3] <= MAX_GAP for i in range(len(sorted_combo)-1)):
                    prox_combos.append(combo)
                    prox = True
            stranded = False
            strand_combos = []
            for combo in prox_combos: 
                strands = [h[4] for h in combo]
                if len(set(strands)) == 1: # if all strands are the same: 
                    strand_combos.append(combo)
                    stranded = True
            order = False
            order_combos = [] 
            for combo in strand_combos: 
                strand = combo[0][4] # strand for this combo
                wrong_order = False
                if strand == '+': 
                    for i in range(0, len(combo)-1): 
                        elem = combo[i]
                        next_elem = combo[i+1]
                        if elem[2] > next_elem[2]:  
                            wrong_order = True
                            break
                elif strand == '-':
                    for i in range(0, len(combo)-1): 
                        elem = combo[i]
                        next_elem = combo[i+1]
                        if elem[2] < next_elem[2]:  
                            wrong_order = True
                            break
                if wrong_order == False:
                    order_combos.append(combo) 
                    order = True
            if len(order_combos) > 0:
                successful_combos.append(order_combos)
            else: 
                if prox == True: 
                    orf_prox = True
                if stranded == True: 
                    orf_strand = True
        if len(successful_combos) > 0:
            print(orf, "INTACT", successful_combos) 
            orf_conserved_anywhere[orf] = True
        else: 
            print(orf, 'failed: proximity:', orf_prox, 'strand:', orf_strand, 'order', orf_order)

for orf in orf_focal_ranges:
    if orf_conserved_anywhere[orf]:
        continue # write nothing
    readid = orf_to_read_dict[orf]
    locusrowidx = reads_to_locus_dict[readid]
    prot_ranges = sorted(orf_focal_ranges[orf], key=lambda r: r[0])
    prot_ids = sorted(extract_name(r[2]) for r in prot_ranges)
    orf_protids = ','.join(prot_ids)
    locus_dict[locusrowidx].append([readid, orf_protids]) # go from orf to read id, readid to locus, append readid to locus as supporting a chimeric orf


with open(reads_to_locus) as f:
    read_to_locus_rows = [line.rstrip('\n').split('\t') for line in f]

outfile = open(sys.argv[3] + '_chimeric_byprot_bylocus_DAFFSPECIFIC', 'w')

for locus in locus_dict: 
    locusrow = read_to_locus_rows[locus] # value is idx
    prefix = locusrow[:-1]
    readconfigs = locus_dict[locus]
    configs = {} # config:readcount
    for config in readconfigs: 
        readconfig = config[1]
        if readconfig not in configs: 
            configs[readconfig] = 1
        else: 
            configs[readconfig] += 1
    for config in configs:
        outfile.write('\t'.join(prefix) + '\t' + 'ORF_CHIMERA' + '\t' + config + '\t' + str(configs[config]) + '\n') 
outfile.close() 


