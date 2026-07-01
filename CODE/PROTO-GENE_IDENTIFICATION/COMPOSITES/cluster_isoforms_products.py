# prints one line for each locus/putative orf combination 
# that is, integrates over isoforms to the extent that they encode the same general protein
# but doesn't to the extent that they encode different proteins or combinations

import sys
import numpy as np 
import re 

#infile = 'allfeature_reads_F_out_gffinfo_TLREADSONLY_analysis_withreads_CODINGORFS_DAFF_SPECIFIC'
infile = sys.argv[1]

with open(infile) as f:
    isolist = [line.rstrip('\n').split('\t') for line in f if line.strip()]

locus_lists = []  # array of arrays; each array contains superset of each fragment in reads that have fragments in common
locus_list_readsperorf = []  # ith array corresponds to ith locus list; components are [type, orf structure, readcount]

for row in isolist:
    foundlocusidx = None
    covinfo = []
    fraglist = []

    for i in range(len(row)):
        entry = row[i]
        frag_re = re.compile(r'_\d+_\d+$')
        if 'melanogaster' in entry and frag_re.search(entry):  # marks a fragment ID with unique coordinates
            fraglist.append(entry)
            if foundlocusidx is None:
                for j in range(len(locus_lists)):
                    if entry in locus_lists[j]:
                        foundlocusidx = j
                        break
        if 'ORF_' in entry:
            covinfo.append([row[i], row[i+1], int(row[i+2])])

    if foundlocusidx is not None:
        samelocus = locus_lists[foundlocusidx]
        for entry in fraglist:
            if entry not in samelocus:
                samelocus.append(entry)
        readsperorf = locus_list_readsperorf[foundlocusidx]
        for currorf in covinfo:
            found = False
            for orf in readsperorf:
                if currorf[0] == orf[0] and currorf[1] == orf[1]:
                    orf[2] += currorf[2]
                    found = True
                    break
            if not found:
                readsperorf.append(currorf)
    else:
        locus_lists.append(fraglist)
        locus_list_readsperorf.append(covinfo)


def parse_orf_shortnames(orfinfo):  # takes names of proteins in this format that are in the orf
    # example: argonaute-3[Drosophila-melanogaster],malate-dehydrogenase-1[Drosophila-melanogaster] 
    prots = re.split(',', orfinfo)
    names = []
    for prot in prots:
        if '-[' in prot:
            name = re.split(r'-\[', prot)[0]
        else: 
            name = re.split(r'\[', prot)[0]
        names.append(name.replace('_', '-'))
    return(names)

def parse_frag_shortnames(fraginfo): # takes names of proteins in fragments in the transcript 
    if '-[' in fraginfo:
        name = re.split(r'-\[', fraginfo)[0]
    else: 
        name = re.split(r'\[', fraginfo)[0]
    return(name.replace('_', '-'))

orftypes = ['ORF_SINGLETON', 'ORF_CHIMERA', 'ORF_NEW_EXON']
orfreadnums = [[] for _ in range(len(orftypes))]  # by type; array of ararys
orfotherfrags = [[] for _ in range(len(orftypes))] # by type; array of arrays; 0 for no other fragments on tx, 1 for other fragments on tx
with open(infile + '_iso_collapsed', 'w') as outfile:
    for i in range(len(locus_lists)):
        fragments = locus_lists[i]
        covinfo = locus_list_readsperorf[i]
        row_parts = list(fragments)
        locus_other_frags = False # to avoid double-counting other frags (if encodes two proteins and one has another frag, the other does too)
        for orf in covinfo:
            other_frags = False
            orf_prots = parse_orf_shortnames(orf[1]) # list of proteins in the orf
            for frag in row_parts:
                 if parse_frag_shortnames(frag) not in orf_prots: 
                     other_frags = True
            if other_frags == True and locus_other_frags == False:
                #print(orf[1])
                locus_other_frags = True
                orfotherfrags[orftypes.index(orf[0])].append(1)
            if other_frags == False: 
                orfotherfrags[orftypes.index(orf[0])].append(0)
            orfinfo = []
            orfinfo.append(orf[0])
            orfinfo.append(orf[1])
            orfinfo.append(str(orf[2]))
            orfreadnums[orftypes.index(orf[0])].append(orf[2])
            outfile.write('\t'.join(orfinfo) + '\t')
            outfile.write('\t'.join(row_parts) + '\n')

            #row_parts.append(orf[0])           # ORF type (e.g. ORF_SINGLETON)
            #row_parts.append(orf[1])           # ORF structure / gene name
            #row_parts.append(str(orf[2]))      # read count

print(infile + '_iso_collapsed')
print('total loci: ', len(locus_lists))
for i in range(0, len(orftypes)): 
    print(orftypes[i], ' number: ', len(orfreadnums[i]))
    print(' average reads: ', np.mean(orfreadnums[i]))
    print(' median reads: ', np.median(orfreadnums[i]))
    print(' fraction with another fragment: ', np.sum(orfotherfrags[i])/len(orfotherfrags[i]))


