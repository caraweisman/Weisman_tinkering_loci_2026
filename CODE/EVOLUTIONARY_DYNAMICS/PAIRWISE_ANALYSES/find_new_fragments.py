import numpy as np 
import sys
import re

#usage: python (code.py) (blast output, self only) (focal mapping file, esl) (conversion file, esl) 
# (fragment gff, TL, focal) (fragment gff, target) (species)

# input files: 
# 1. blast output file between corresponding TL segments in focal and target
# 2. focal mapping file to map the SUBregion coordinates in affinis to their genomic coordinates
# 2. conversion file to map the region corresponding to the focal species to the true coordinates in the target species
# 3. gff of fragments from focal species (only TL in this case) 
# 4. gff of fragments from target species (all) 

blastout = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t') 
focal_mapping = np.genfromtxt(sys.argv[2], dtype=str, delimiter='\t') # this sucks but is to convert the subregion in affinis
# which are numbered by the original overall TL back to their coordinates in the genome
conversion = np.genfromtxt(sys.argv[3], dtype=str, delimiter='\t') # this sucks but is the mapping between the regions in the
# target fasta which are NAMED FOR THE FOCAL REGION and their actual coordinates in the focal genome. 
# it's the sfetch file for the target.
focal_gff = np.genfromtxt(sys.argv[4], dtype=str, delimiter='\t') 
target_gff = np.genfromtxt(sys.argv[5], dtype=str, delimiter='\t') 
species = sys.argv[6]
self_blastout=np.genfromtxt('Daff_reconstruction_seqs.fna_vs_Daff_reconstruction_seqs.fna_blastout_selfonly_nodiag', dtype=str, delimiter='\t')
target_self_blastout = np.genfromtxt(species + '_reconstruction_seqs.fna_vs_' + species + '_reconstruction_seqs.fna_blastout_selfonly_nodiag', dtype=str, delimiter='\t')

def extract_name_short(desc):
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

def extract_name(desc):
    ID = re.split(';', desc)[0]
    body = re.split('=', ID, maxsplit=1)[1]
    name = re.split('-\\[', body)[0]
    return name



# make dictionary of focal-target region conversion to get target contigs and coordinates
conversion_dict = {}
for row in conversion:
    locid = row[0]
    start = int(row[1])
    end = int(row[2])
    contig = row[3]
    conversion_dict[locid] = [contig, start, end]

# make dictionary of mapping from focal SUBregions to their actual genomic coordinates
focal_mapping_dict = {}
for row in focal_mapping: 
    locid = row[0]
    start = int(row[1])
    end = int(row[2])
    contig = row[3]
    focal_mapping_dict[locid] = [contig, start, end]


gff_to_blast_strand_dict = {'+':'plus', '-':'minus'} # different encodings; gff has +/-, blast has plus/mins

# make dictionary of fragments by contig
# sort first by ocntig and then by position - so that relative order can be assessed 
focal_gff = focal_gff[np.lexsort((focal_gff[:,3].astype(int), focal_gff[:,0]))]
target_gff = target_gff[np.lexsort((target_gff[:,3].astype(int), target_gff[:,0]))]
# first for focal species
focal_frag_dict = {}
for frag in focal_gff: 
    start = int(frag[3])
    end = int(frag[4])
    strand = gff_to_blast_strand_dict[frag[6]]
    contig = frag[0]
    name = extract_name(frag[8])
    if contig not in focal_frag_dict: 
        focal_frag_dict[contig] = []
    focal_frag_dict[contig].append([start, end, name, strand])
# then for target species
target_frag_dict = {}
for frag in target_gff: 
    start = int(frag[3])
    end = int(frag[4])
    strand = gff_to_blast_strand_dict[frag[6]]
    contig = frag[0]
    name = extract_name(frag[8])
    if contig not in target_frag_dict: 
        target_frag_dict[contig] = []
    target_frag_dict[contig].append([start, end, name, strand])


# dictionary that will output how many distinct regions in the focal / target with fragments there are in the target/focal
# because the rest of the script only judges flux based on subregion
num_subregions_with_frags_dict = {}
for row in conversion: # = subregion
    rawregion = row[0]
    fullregionid = rawregion.rsplit('-', 1)[0] # the whole region, not the subregion
    if fullregionid not in num_subregions_with_frags_dict: # don't have the raw region list here; have the file but don't want to have it as an input
        num_subregions_with_frags_dict[fullregionid] = [[],[]] # two lists
        # first list: number of fragments on each subregion in the focal (by definition, noncontiguous with subregion in target) 
        # second list: number of fragments on each subregion in target 
        # each looks like eg 1,1,2; three subregions with 1, 1, and 2 fragments respectively

# make dictionary of all unique regions in the blast output file 
# and what fragments are in them
# do for focal and target at the same time (because they're both in the conversion file)
#### put it in relative coordinates: measured from the start of the region
# so that it can be easily mapped to the blast file
focal_regions_fragments = {}
target_regions_fragments = {}
for row in conversion: # = subregion
    # first focal fragments
    # coordinates for region
    focalid = row[0]
    numfocalfrags = 0
    fullregionid = focalid.rsplit('-', 1)[0] # the whole region, not the subregion
    focalcontig, focalstart,focalstop =  focal_mapping_dict[focalid]
    for frag in focal_frag_dict.get(focalcontig, []): # I didn't think this would be necessary, but the HMM has segmented some contigs as TLs with no fragments.........
        # these should be in ascending order of position
        # coordinates for fragment
        fragstart = frag[0]
        fragstop = frag[1]
        fragname = frag[2]
        strand = frag[3]
        fraglen = fragstop - fragstart
        overlap = max(0, min(focalstop, fragstop) - max(focalstart, fragstart))
        if overlap >= float(fraglen)/2: 
            if focalid not in focal_regions_fragments: 
                focal_regions_fragments[focalid] = []
            # append relative coordinates if they don't overlap with something else
            # this only happens for short fragments below the overlap threshold in the annotation pipeline
            relstart = fragstart - focalstart
            relstop = fragstop - focalstart
            fragoverlap = False
            for existing in focal_regions_fragments[focalid]:
                existingstart = existing[0]
                existingstop = existing[1]
                existinglen = existingstop - existingstart
                existingoverlap = max(0, min(relstop, existingstop) - max(relstart, existingstart))
                if existingoverlap >= float(min(fraglen, existinglen))/2:
                    fragoverlap = True
                    break
            if fragoverlap == False:
                focal_regions_fragments[focalid].append([relstart, relstop, fragname, strand])
                numfocalfrags += 1
    num_subregions_with_frags_dict[fullregionid][0].append(numfocalfrags) #0th entry is focal vs target
    # now target fragments; should also be in ascending order of position
    targetcontig, targetstart, targetstop = conversion_dict[focalid] # use conversion dict which was already made
    numtargetfrags = 0
    for frag in target_frag_dict.get(targetcontig, []):
        fragstart = frag[0]
        fragstop = frag[1]
        fragname = frag[2]
        strand = frag[3]
        fraglen = fragstop - fragstart
        overlap = max(0, min(targetstop, fragstop) - max(targetstart, fragstart))
        if overlap >= float(fraglen)/2: 
            if focalid not in target_regions_fragments: 
                target_regions_fragments[focalid] = []
            # append relative coordinates
            relstart = fragstart - targetstart
            relstop = fragstop - targetstart
            fragoverlap = False
            for existing in target_regions_fragments[focalid]:
                existingstart = existing[0]
                existingstop = existing[1]
                existinglen = existingstop - existingstart
                existingoverlap = max(0, min(relstop, existingstop) - max(relstart, existingstart))
                if existingoverlap >= float(min(fraglen, existinglen))/2:
                    fragoverlap = True
                    break
            if fragoverlap == False:
                target_regions_fragments[focalid].append([relstart, relstop, fragname, strand]) 
                numtargetfrags += 1
    num_subregions_with_frags_dict[fullregionid][1].append(numtargetfrags) #1st entry is target vs focal

# write outfile for number of subregions with fragments 
outfile0 = open ('discontiguous_region_fragments', 'w')
for region in num_subregions_with_frags_dict: 
    outfile0.write(region + '+') 
    focal_subregion_counts = num_subregions_with_frags_dict[region][0]
    for subregion in focal_subregion_counts: 
        outfile0.write('\t' + str(subregion))
    outfile0.write('\n') 
    outfile0.write(region + '-') 
    target_subregion_counts = num_subregions_with_frags_dict[region][1]
    for subregion in target_subregion_counts: 
        outfile0.write('\t' + str(subregion))
    outfile0.write('\n') 
outfile0.close()


# make dictionary of blast results for each region so don't have to search a gazillion times
region_blastouts = {}
for row in blastout: 
    regionid = row[0]
    if regionid not in region_blastouts: 
        region_blastouts[regionid] = []
    region_blastouts[regionid].append(row)

# also make one for focal vs focal search for in-duplicates
self_region_blastouts = {}
for row in self_blastout: 
    regionid = row[0]
    if regionid not in self_region_blastouts: 
        self_region_blastouts[regionid] = []
    self_region_blastouts[regionid].append(row)

# and one for target vs target search for in-duplicates
target_self_region_blastouts = {}
for row in target_self_blastout: 
    regionid = row[0]
    if regionid not in target_self_region_blastouts: 
        target_self_region_blastouts[regionid] = []
    target_self_region_blastouts[regionid].append(row)



# create dictionary of all TL subregions and what fragments are in them
missing_fragment_dict = {}
missing_duplicate_dict = {} # to print duplicates separately; previous version of code had no way to distinguish

outfile = open('fragment_arrangement_' + species, 'w')
outfile.write('# top row for each region is fragments in focal rearranged in target; \n')
outfile.write('# bottom row for each region is fragments in target rearranged in focal; \n')
for region in conversion:
    # first check fragments in query to see if they're covered in the regions in the blast alignment to target
    focalid = region[0]
    # initialize entry in dict
    # format: key is SUBregion ID (can stitch back together at the end),
    # value is two arrays, the first with fragments from the focal species
    # and the second with fragments from the target species
    missing_fragment_dict[focalid] = [[],[]]
    missing_duplicate_dict[focalid] = [[],[]]
    blasthits = region_blastouts.get(focalid, []) # in case no blast hits for this region 
    exp_focal_fragments = focal_regions_fragments.get(focalid, []) # fragments expected in this region from gff
    exp_target_fragments = target_regions_fragments.get(focalid, []) # fragments expected in this region from gff
    expfocalfragorder = [] # array: (fragid, strand) 2ple. 
    targetfocalfragorder = [] # array (start, fragid, strand) 4ple; tracks realized order of focal fragments IN TARGET. use start only - i think this works?
    # make array to keep track of which target coordinates (normalized to region) have been 'used' by each hit, to detect duplications in focal 
    used_target_coords = []
    # make array to keep track of which focal coordinates (norm to region) have been used by each hit counted as a GAIN in the target, to keep track of which hits are in-duplicates
    used_focal_coords_gainfrags = []
    for fragment in exp_focal_fragments: 
        fragstart = fragment[0]
        fragend = fragment[1]
        fragname = fragment[2]
        fragstrand = fragment[3]
        fraglen = fragend-fragstart
        found = False
        blast_passed = False
        # fragments are sorted in order 
        # exclude fragments from any same gene in the other region; conservative 
        #if fragment in [f[2] for f in exp_target_fragments]:
            #found = True
        for z in range(0, len(blasthits)):
            hit = blasthits[z]  
            # query start/end in blast output
            hitstart = min(int(hit[2]), int(hit[3]))  
            hitend = max(int(hit[2]), int(hit[3]))
            targethitstart = min(int(hit[6]), int(hit[7])) 
            hitstrand = hit[11] 
            if hitstrand == 'plus':
                targetfragstart = int(hit[6]) + (fragstart - hitstart) # project start/end coords of frag onto target; not precise, but probably close
                targetfragend = int(hit[6]) + (fragend - hitstart)
            else: 
                targetfragstart = int(hit[6]) - (fragstart - hitstart)
                targetfragend = int(hit[6]) - (fragend - hitstart)
            if hitstrand == 'plus':
                targetfragstrand = fragstrand
            else: 
                targetfragstrand = 'minus' if fragstrand == 'plus' else 'plus'
            overlap = max(0, min(fragend, hitend) - max(fragstart, hitstart))
            if targetfragstart > targetfragend:
                targetfragstart, targetfragend = targetfragend, targetfragstart
            # clamp projection to hit's target boundaries to prevent overlap from bad blast inference due to hsps ending and inferred coords running over 
            target_hit_lo = min(int(hit[6]), int(hit[7]))
            target_hit_hi = max(int(hit[6]), int(hit[7]))
            targetfragstart = max(targetfragstart, target_hit_lo)
            targetfragend = min(targetfragend, target_hit_hi)
            puttargetfraglen = targetfragend-targetfragstart # use clamped length
            if overlap >= float(fraglen)*0.25: 
                blast_passed = True
                if not any(max(0, min(targetfragend, e) - max(targetfragstart, s)) >= puttargetfraglen*0.5 for s, e in used_target_coords):
                    found = True 
                    # frags are listed in order in the gff and thus the dict; easy
                    # append information about position and strand from both target and focal only if found in target
                    expfocalfragorder.append([fragname, fragstrand])
                    targetfocalfragorder.append([targetfragstart, fragname, targetfragstrand])
                    used_target_coords.append([targetfragstart, targetfragend])
                    break # takes best nucleotide-level blast hit between genomes that finds each fragment
        if found == False:
            if blast_passed == False:
                selfhits = False
                already_gained = False
                selfblasthits = self_region_blastouts.get(focalid, [])
                for selfhit in selfblasthits:
                    qstart = min(int(selfhit[2]), int(selfhit[3]))
                    qend = max(int(selfhit[2]), int(selfhit[3]))
                    sstart_raw = int(selfhit[6])
                    send_raw = int(selfhit[7])
                    sub_lo = min(sstart_raw, send_raw)
                    sub_hi = max(sstart_raw, send_raw)
                    selfhitstrand = selfhit[11]
                    overlap = max(0, min(fragend, qend) - max(fragstart, qstart))
                    if overlap >= float(fraglen)*0.25:
                        selfhits = True
                        for usedstart, usedend in used_focal_coords_gainfrags:
                            used_in_subject = max(0, min(usedend, sub_hi) - max(usedstart, sub_lo))
                            used_len = usedend - usedstart
                            if used_in_subject < used_len * 0.5:
                                continue
                            if selfhitstrand == 'plus':
                                proj_start = qstart + (usedstart - sub_lo)
                                proj_end = qstart + (usedend - sub_lo)
                            else:
                                proj_start = qstart + (sub_hi - usedend)
                                proj_end = qstart + (sub_hi - usedstart)
                            if proj_start > proj_end:
                                proj_start, proj_end = proj_end, proj_start
                            proj_start = max(proj_start, qstart)
                            proj_end = min(proj_end, qend)
                            proj_len = proj_end - proj_start
                            if proj_len <= 0:
                                continue
                            frag_proj_overlap = max(0, min(fragend, proj_end) - max(fragstart, proj_start))
                            if frag_proj_overlap >= proj_len * 0.5:
                                already_gained = True
                                break
                        if already_gained:
                            break
                if selfhits:
                    missing_duplicate_dict[focalid][0].append(fragname)
                    if not already_gained:
                        missing_fragment_dict[focalid][0].append(fragname)
                        used_focal_coords_gainfrags.append([fragstart, fragend])
                else:
                    missing_fragment_dict[focalid][0].append(fragname)
                    used_focal_coords_gainfrags.append([fragstart, fragend])
            elif blast_passed == True:
                missing_duplicate_dict[focalid][0].append(fragname)
    # check if all fragments in region are in the same relative order and orientation 
    # sort order in fragments
    targetfocalfragorder.sort(key=lambda x: int(x[0]))
    print(focalid)
    print(expfocalfragorder)
    print(targetfocalfragorder)
    # remove start position now that it's sorted
    targetfocalfragorder = [x[1:] for x in targetfocalfragorder]
    consconfig = False
    consconfigplus = True # only looks for rearrangements; ignores insertions/deletions, which earlier branch detects
    consconfigminus = True
    for i in range(0, len(expfocalfragorder)):
        if expfocalfragorder[i][0] == targetfocalfragorder[i][0] and expfocalfragorder[i][1] == targetfocalfragorder[i][1]: 
            continue
        else: 
            consconfigplus = False
    for i in range(0, len(expfocalfragorder)):
        if expfocalfragorder[i][0] == targetfocalfragorder[len(targetfocalfragorder) - 1 - i][0] and expfocalfragorder[i][1] != targetfocalfragorder[len(targetfocalfragorder) - 1 - i][1]: 
            continue
        else:
            consconfigminus = False
    if consconfigminus == True or consconfigplus == True: 
        consconfig = True
    outfile.write(focalid + '\t')
    if consconfig == False: 
        outfile.write('REARRANGED' + '\n')
    if consconfig == True:
        outfile.write('CONSERVED' + '\n')
    # now check fragments in target
    exptargetfragorder = [] # array: (fragid, strand) 2ple. 
    focaltargetfragorder = [] # array (start, fragid, strand) 3ple; tracks realized order of target fragments IN FOCAL. use start only - i think this works?
    used_focal_coords = [] # detect which inferred regions in the focal region have been 'used' to detect duplications
    used_target_coords_gainfrags = [] # track which in-duplication(s) should also be called as a gain
    for fragment in exp_target_fragments: 
        fragstart = fragment[0]
        fragend = fragment[1]
        fragname = fragment[2]
        fragstrand = fragment[3] # strand in TARGET
        fraglen = fragend-fragstart
        found = False
        blast_passed = False
        # exclude fragments from any same gene in the other region; conservative 
        #if fragname in [f[2] for f in exp_focal_fragments]:
        #    found = True
        for z in range(0, len(blasthits)): 
            hit = blasthits[z]
            # target start/end in blast output
            hitstart = min(int(hit[6]), int(hit[7])) 
            hitend = max(int(hit[6]), int(hit[7]))
            hitstrand = hit[11] 
            overlap = max(0, min(fragend, hitend) - max(fragstart, hitstart))
            if hitstrand == 'plus':
                focalfragstart = int(hit[2]) + (fragstart - hitstart)
                focalfragend = int(hit[2]) + (fragend - hitstart)
            else: 
                focalfragstart = int(hit[3]) - (fragstart - hitstart)
                focalfragend = int(hit[3]) - (fragend - hitstart)
            if hitstrand == 'plus':
                targetfragstrand = fragstrand
            else: 
                targetfragstrand = 'minus' if fragstrand == 'plus' else 'plus'
            if focalfragstart > focalfragend:
                focalfragstart, focalfragend = focalfragend, focalfragstart
            focal_hit_lo = min(int(hit[2]), int(hit[3]))
            focal_hit_hi = max(int(hit[2]), int(hit[3]))
            focalfragstart = max(focalfragstart, focal_hit_lo)
            focalfragend = min(focalfragend, focal_hit_hi)
            putfocalfraglen = focalfragend - focalfragstart
            if overlap >= float(fraglen)*0.25: 
                blast_passed = True
                if not any(max(0, min(focalfragend, e) - max(focalfragstart, s)) >= putfocalfraglen*0.5 for s, e in used_focal_coords):
                    found = True
                    exptargetfragorder.append([fragname, fragstrand])
                    focaltargetfragorder.append([focalfragstart, fragname, targetfragstrand])
                    used_focal_coords.append([focalfragstart, focalfragend])
                    break # takes best nucleotide-level blast hit between genomes that finds each fragment
        if found == False:
            if blast_passed == False:
                selfhits = False
                already_gained = False
                selfblasthits = target_self_region_blastouts.get(focalid, [])
                for selfhit in selfblasthits:
                    qstart = min(int(selfhit[2]), int(selfhit[3]))
                    qend = max(int(selfhit[2]), int(selfhit[3]))
                    sstart_raw = int(selfhit[6])
                    send_raw = int(selfhit[7])
                    sub_lo = min(sstart_raw, send_raw)
                    sub_hi = max(sstart_raw, send_raw)
                    selfhitstrand = selfhit[11]
                    overlap = max(0, min(fragend, qend) - max(fragstart, qstart))
                    if overlap >= float(fraglen)*0.25:
                        selfhits = True
                        for usedstart, usedend in used_target_coords_gainfrags:
                            used_in_subject = max(0, min(usedend, sub_hi) - max(usedstart, sub_lo))
                            used_len = usedend - usedstart
                            if used_in_subject < used_len * 0.5:
                                continue
                            if selfhitstrand == 'plus':
                                proj_start = qstart + (usedstart - sub_lo)
                                proj_end = qstart + (usedend - sub_lo)
                            else:
                                proj_start = qstart + (sub_hi - usedend)
                                proj_end = qstart + (sub_hi - usedstart)
                            if proj_start > proj_end:
                                proj_start, proj_end = proj_end, proj_start
                            proj_start = max(proj_start, qstart)
                            proj_end = min(proj_end, qend)
                            proj_len = proj_end - proj_start
                            if proj_len <= 0:
                                continue
                            frag_proj_overlap = max(0, min(fragend, proj_end) - max(fragstart, proj_start))
                            if frag_proj_overlap >= proj_len * 0.5:
                                already_gained = True
                                break
                        if already_gained:
                            break
                if selfhits:
                    missing_duplicate_dict[focalid][1].append(fragname)
                    if not already_gained:
                        missing_fragment_dict[focalid][1].append(fragname)
                        used_target_coords_gainfrags.append([fragstart, fragend])
                else:
                    missing_fragment_dict[focalid][1].append(fragname)
                    used_target_coords_gainfrags.append([fragstart, fragend])
            elif blast_passed == True:
                missing_duplicate_dict[focalid][1].append(fragname)

    # check if all fragments in region are in the same relative order and orientation 
    # sort order in fragments
    focaltargetfragorder.sort(key=lambda x: int(x[0]))
    print(focalid)
    print(exptargetfragorder)
    print(focaltargetfragorder)
    # remove start position now that it's sorted
    focaltargetfragorder = [x[1:] for x in focaltargetfragorder]
    consconfig = False
    consconfigplus = True # only looks for rearrangements; ignores insertions/deletions, which earlier branch detects
    consconfigminus = True
    for i in range(0, len(exptargetfragorder)):
        if exptargetfragorder[i][0] == focaltargetfragorder[i][0] and exptargetfragorder[i][1] == focaltargetfragorder[i][1]: 
            continue
        else: 
            consconfigplus = False
    for i in range(0, len(exptargetfragorder)):
        if exptargetfragorder[i][0] == focaltargetfragorder[len(focaltargetfragorder) - 1 - i][0] and exptargetfragorder[i][1] != focaltargetfragorder[len(focaltargetfragorder) - 1 - i][1]: 
            continue
        else: 
            consconfigminus = False
    if consconfigminus == True or consconfigplus == True: 
        consconfig = True
    outfile.write(focalid + '\t')  
    if consconfig == False: 
        outfile.write('REARRANGED' + '\n')
    if consconfig == True:
        outfile.write('CONSERVED' + '\n')
outfile.close()
            
outfile = open('fragment_flux_' + species, 'w') 
outfile.write('# top row for each region is fragments in focal missing in target; \n')
outfile.write('# bottom row for each region is fragments in target missing in focal; \n')
for region in missing_fragment_dict: 
    missing_frags_focal = missing_fragment_dict[region][0]
    missing_frags_target = missing_fragment_dict[region][1]
    outfile.write(region)
    for frag in missing_frags_focal: 
        outfile.write('\t' + frag) 
    outfile.write('\n')
    outfile.write(region)
    for frag in missing_frags_target: 
        outfile.write('\t' + frag) 
    outfile.write('\n')    
outfile.close() 
    

outfile = open('duplicate_flux_' + species, 'w') 
outfile.write('# top row for each region is fragments in focal missing in target; \n')
outfile.write('# bottom row for each region is fragments in target missing in focal; \n')
for region in missing_duplicate_dict: 
    missing_frags_focal = missing_duplicate_dict[region][0]
    missing_frags_target = missing_duplicate_dict[region][1]
    outfile.write(region)
    for frag in missing_frags_focal: 
        outfile.write('\t' + frag) 
    outfile.write('\n')
    outfile.write(region)
    for frag in missing_frags_target: 
        outfile.write('\t' + frag) 
    outfile.write('\n')    
outfile.close() 
    













