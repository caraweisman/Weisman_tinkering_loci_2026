import numpy as np
import sys
from collections import Counter
import re 

#speciestree = [['Dath', 'Dalg'], 'Dazt', 'Dhel']
#species = ['Dath', 'Dalg', 'Dazt', 'Dhel']
species = ['Dhel', 'Dazt', 'Dalg', 'Dath']
#speciestree = [[0,1],2,3]
speciestree = [0,1,2,3]


def long_to_short_name(name): 
	desc = '-'.join(re.split('-', name)[3:])
	isofree = re.split(',', desc)[0]
	return isofree


regions = np.genfromtxt('TL_sfetch', dtype=str, delimiter='\t', usecols=(0), comments='#')

resultfiles = ['duplicate_flux_' + s for s in species]

# initialize dictionary that will hold compiled results from all species
region_flux_dict = {}
for region in regions: 
	region_flux_dict[region] = [] # for every region, an entry in the dict
	for s in species: 
		region_flux_dict[region].append([[],[]]) # for every region, the dict entry gets one array per species, consisting of two subarrays: things present in focal but missing
# in this species, and present in this species and missin from focal


# check how many subregions WITH FRAGS are in each region
# and populate dictionary for missing frags
speciesidx = 0 # increment this!!!!!!!!! 
for file in resultfiles: 
	seenregions = []
	with open(file) as f:
		for line in f:
			row = line.rstrip('\n').split('\t')
			rawregion = row[0]
			region = rawregion.rsplit('-', 1)[0]  # to remove suffix -1, -2, etc for partitioning TLs into alignable regions
			hits = row[1:]
			if rawregion not in seenregions: # add at the end of this!!!!!!!  
				seenregions.append(rawregion)
				for hit in hits: 
					shortname = long_to_short_name(hit)
					region_flux_dict[region][speciesidx][0].append(hit) 
			else: 
				for hit in hits: 
					shortname = long_to_short_name(hit)
					region_flux_dict[region][speciesidx][1].append(hit) 
	speciesidx += 1


# analyze dictionary

# make new dictionary of 'compiled' fluxes, taking phylogeny into account
# per region, six entries: one "plus/minus" per taxon: in focal and not in target, and vice versa
# contains names of fragments, which will be enforced to be unique
phylo_fluxes = {} 
for region in regions: 
	phylo_fluxes[region] = []
	for taxon in speciestree: 
		phylo_fluxes[region].append([[],[]]) # plus, minus

# this is a little convoluted, but the idea here is to give an accounting of the total flux at the locus over the group. 
# each new fragment at the locus gets counted ONCE, in some species. so rough interpretation is: variation in fragment content across the group.
# a mild phylogenetic interpretation based on which species it's listed as being in.
# for fragments present in focal affinis and absent in one or more outgroup species, it gets counted in the earliest-diverging species in which it is absent. 
# for fragments present in an outgroup and absent in affinis, it gets counted in the earliest diverging species in which it is present. 
# second interpretation is easier: age of the fragment; others explainable by losses. 
# first interpretation depends more on details. 

for region in regions: 
	fluxes = region_flux_dict[region]
	# integrate over all species in taxon
	for taxon in speciestree: # taxa in speciestree are numbers whose positions correspond to the positions of species in species list and therefore in flux dict subarrays
		taxonidx = speciestree.index(taxon)
		members = taxon if isinstance(taxon, list) else [taxon] 
		focal_frags = [] # in focal but not target
		target_frags = [] # in target but not focal
		for s in members: 
			curr_focal_frags = fluxes[s][0]
			for frag in curr_focal_frags: 
				focal_frags.append(frag) 
			curr_target_frags = fluxes[s][1]
			for frag in curr_target_frags: 
				target_frags.append(frag) 
		previous_indices = list(range(taxonidx))
		old_focal_frags = []
		old_target_frags = []
		for prevtaxa in previous_indices: 
			old_focal_frags += phylo_fluxes[region][prevtaxa][0]
			old_target_frags += phylo_fluxes[region][prevtaxa][1]
		focal_counts = Counter(focal_frags)
		old_focal_counts = Counter(old_focal_frags)
		for fragname, count in focal_counts.items():
			new_count = count - old_focal_counts.get(fragname, 0) # how many more fragments are in the focal vs current target than were previously seen in focal vs any previous target
			for _ in range(max(0, new_count)): # if current target has fewer relative to focal than previous, add nothing (0 max)
				phylo_fluxes[region][taxonidx][0].append(fragname) # otherwise, add one entry for each new copy 
		target_counts = Counter(target_frags)
		old_target_counts = Counter(old_target_frags)
		for fragname, count in target_counts.items():
			new_count = count - old_target_counts.get(fragname, 0)
			for _ in range(max(0, new_count)):
				phylo_fluxes[region][taxonidx][1].append(fragname)


#outfile = open('TL_dupfrags_flux_numbers', 'w')
#outfile.write('#Region \t Dhel \t Dazt \t Dalg \t Dath \t Total \n')

outfile2 = open('Total_duplicatedfrags_flux_numbers', 'w') # adds focal/target
outfile2.write('#Region \t Total \n')

for region in phylo_fluxes:
	totalnumfocalfrags = 0
	totalnumtargetfrags = 0 
	#outfile.write(region + '+')
	for taxon in phylo_fluxes[region]: 
		focals = taxon[0]
		#fraglist = ','.join(focals) 
		#outfile.write('\t' + fraglist) 
		#outfile.write('\t' + str(len(focals)))
		totalnumfocalfrags += len(focals)
	#outfile.write('\t' + str(totalnumfocalfrags))
	#outfile.write('\n')
	#outfile.write(region + '-') 
	for taxon in phylo_fluxes[region]: 
		targets = taxon[1] 
		#fraglist = ','.join(targets) 
		#outfile.write('\t' + fraglist) 
		#outfile.write('\t' + str(len(targets)))
		totalnumtargetfrags += len(targets)
	#outfile.write('\t' + str(totalnumtargetfrags))
	#outfile.write('\n')
	totalfrags = totalnumfocalfrags + totalnumtargetfrags
	outfile2.write(region + '\t' + str(totalfrags) + '\n')

#outfile.close()
outfile2.close()


outfile = open('Dupfrags_flux_IDs', 'w')
outfile.write('#Region \t Dhel \t Dazt \t Dalg \t Dath \n')
for region in phylo_fluxes: 
	outfile.write(region + '+')
	for taxon in phylo_fluxes[region]: 
		focals = taxon[0]
		fraglist = ','.join(focals) 
		outfile.write('\t' + fraglist) 
		#outfile.write('\t' + str(len(focals)))
	outfile.write('\n')
	outfile.write(region + '-') 
	for taxon in phylo_fluxes[region]: 
		targets = taxon[1] 
		fraglist = ','.join(targets) 
		outfile.write('\t' + fraglist) 
		#outfile.write('\t' + str(len(targets)))
	outfile.write('\n')
outfile.close()


