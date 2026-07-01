import numpy as np
from statsmodels.stats.multitest import multipletests
import sys
import scipy
from scipy.stats import binomtest, binom, poisson
from math import lgamma, log
import re
 
## ***** counts each unique gene at each locus once, regardless of fragment boundaries
# duplications of the same fragment are identified and treated in the FRAG_DUP directory; these are not duplicatinos
# multiple pieces of the same gene, conservatively,  result from different exons of the same gene or fragmentation or loss of parts after gain
# assigns gain of fragment(s) from each unique gene to a branch in the phylogeny
# does a statistical test to see if a given locus has a nonuniform rate of fragment gain across the phylogeny
# updates 5/25: - added chi square test for overall test for heterogeniety; above test is per-branch, this is not. higher power, lower resolution
# - changed from bonferroni correction here to second sweep to FDR correct given the best case p-values as the number of tests
# updated 5/30: changed exclusion of fragments based on coverage to only in the case where their ambiguous state due to coverage actually affects the result; if so, then exclude; if not, then keep
# updated: changed chi square test to eaxct multinomial because of unreliability of p value calculation at low expected counts
# and added a per-locus test for at least one significant branch to help with power 
# updated: locus-level "at least one significant branch" enrichment test now uses exact joint multinomial null instead of independence-assumption union formula

species = ['Dhel', 'Dazt', 'Dalg', 'Dath']
speciestree = [0,1,2,3]

regions = np.genfromtxt('TL_sfetch', dtype=str, delimiter='\t', usecols=(0), comments='#')

resultfiles = ['fragment_flux_' + s for s in species]


def long_to_short_name(name): # flux output files include framgent locations; remove this in conservative assumption that they all stem from the same gain event
        desc = '-'.join(re.split('-', name)[3:])
        isofree = re.split(',', desc)[0]
        return isofree


# helpers for exact multinomial goodness-of-fit test
def enumerate_compositions(n, k): 
	if k == 1: 
		yield (n,)
		return
	for i in range(n+1): 
		for rest in enumerate_compositions(n-i, k-1): 
			yield (i,) + rest

def log_multinomial_pmf(counts, probs, n): 
	logp = lgamma(n+1)
	for c, p in zip(counts, probs): 
		if c > 0: 
			if p <= 0: 
				return -float('inf')
			logp -= lgamma(c+1)
			logp += c * log(p)
	return logp

def chisq_stat(counts, expected): 
	s = 0.0
	for c, e in zip(counts, expected): 
		if e > 0: 
			s += (c - e)**2 / e
	return s

def exact_multinomial_pvalue(observed, probs, n): 
	expected = [n * p for p in probs]
	obs_stat = chisq_stat(observed, expected)
	logps = []
	stats = []
	for comp in enumerate_compositions(n, len(probs)): 
		stats.append(chisq_stat(comp, expected))
		logps.append(log_multinomial_pmf(comp, probs, n))
	logps_arr = np.array(logps)
	stats_arr = np.array(stats)
	mask = stats_arr >= obs_stat - 1e-12
	mx = logps_arr.max()
	pval = np.exp(logps_arr[mask] - mx).sum() / np.exp(logps_arr - mx).sum()
	return float(pval), obs_stat


# initialize dictionary that will hold compiled results from all species
region_flux_dict = {}
for region in regions: 
	region_flux_dict[region] = [] # for every region, an entry in the dict
	for s in species: 
		region_flux_dict[region].append([[],[]]) 
# for every region, the dict entry gets one array per species, consisting of two subarrays: things present in focal but missing
# in this species, and present in this species and missin from focal

# populate dictionary with data from individual species' flux files
speciesidx = 0 
for file in resultfiles: 
	seenregions = [] # hacky way to keep track of the output format from the previous step, in which each locus has two lines, top meaning present in focal/absent in target and bottom vice versa
	with open(file) as f:
		for line in f:
			row = line.rstrip('\n').split('\t')
			rawregion = row[0]
			region = rawregion.rsplit('-', 1)[0]  # to remove suffix -1, -2, etc for partitioning TLs into alignable regions
			hits = row[1:]
			if rawregion not in seenregions: 
				seenregions.append(rawregion)
				for hit in hits: 
					if long_to_short_name(hit) not in region_flux_dict[region][speciesidx][0]:
						region_flux_dict[region][speciesidx][0].append(long_to_short_name(hit)) 
			else: 
				for hit in hits: 
					if long_to_short_name(hit) not in region_flux_dict[region][speciesidx][1]:
						region_flux_dict[region][speciesidx][1].append(long_to_short_name(hit)) 
	speciesidx += 1

# assign branches and lengths
numbranches = 8 # 8 within the phylogeny with lengths; ancestral branch (8) is of indeterminate length
#branchlens = [3.9,1.8,2.1,0.6,0.2,1.3,1.3,1.5] # from BK HOG UCLD tree, PLOS Biology 2024
branchlens = [0.063, 0.041, 0.063, 0.013, 0.005, 0.028, 0.035, 0.046]

# initialize dictionary that will hold fragment assignments per branch per region, to be populated below in loop
region_branch_dict = {}

# load coverage information for each region from minimap orthology search
# file has species in same order as the species here by design (see that code) 
coverage_stats = np.genfromtxt('compiled_coverage', dtype=str, delimiter='\t', skip_header=True)
coverage_dict = {} 
coverage_thresh = 0.95
for row in coverage_stats: 
	region = row[0]
	coverage_vals = row[1:]
	coverage_rulings = []
	for val in coverage_vals: 
		if float(val) < coverage_thresh: 
			coverage_rulings.append(0) 
		else: 
			coverage_rulings.append(1) 
	coverage_dict[region] = coverage_rulings


 
for region in region_flux_dict: 
	coverage_rulings = coverage_dict[region]
	region_branch_dict[region] = [[] for _ in range(numbranches+1)] 
	missing_aff_frag_list = [i[0] for i in region_flux_dict[region]] # this is a list of lists: fragments in each species that are present in affinis and absent in it
	missing_aff_frag_set = {frag for i in region_flux_dict[region] for frag in i[0]} # set: all fragments in affinis missing in at least one species
	for frag in missing_aff_frag_set: 
		absence = [int(frag in s) for s in missing_aff_frag_list] # 1 if absent, 0 if present. only considering gains of unique genes here, so this is appropriate.
		first_present = next((s for s in range(len(coverage_rulings)) if absence[s] == 0 and coverage_rulings[s] == 1), None) # oldest species confidently present
		first_suspect = next((s for s in range(len(coverage_rulings)) if absence[s] == 0 and coverage_rulings[s] == 0), None) # oldest low-cov species marked present (present is the default when region can't be found)
		if first_suspect is not None and (first_present is None or first_suspect < first_present): # only exclude when a suspect species older than the oldest confident presence would change the oldest-presence call
			continue
		oldest = absence.index(0) if 0 in absence else None # finds leftmost (oldest) 0 in array; array is absence, so 0 = present; finds oldest presence.
		# since all frags are definitionally not in all species, assign age based on oldest species; if in affinis and that, assume it was gained on the branch before their divergence.
		if oldest == None:
			branch = 7 # affinis specific
		elif oldest == 0: 
			branch = 8 # in helvetica and affinis; ancestral branch
		elif oldest == 1: 
			branch = 1 # in azteca and affinis; branch leading to dazt-other split
		elif oldest == 2 or oldest == 3:  
			branch = 3 # in alg/ath and affinis; branch leading to ath/alg/aff split
		region_branch_dict[region][branch].append(frag) 
	# repeat the process for fragments present in target species and absent in affinis 
	missing_target_frag_list = [i[1] for i in region_flux_dict[region]] # this is a list of lists: fragments in each species that are absent in affinis and present in a target
	# for each target species, consider all of its fragments missing in affinis
	for i in range(0, len(missing_target_frag_list)): 
		spec_frags = missing_target_frag_list[i]
		for frag in spec_frags: # for each fragment in this species
			presence = [] # this time, if a frag is in another species' list, it means it's present in both (and absent in affinis by construction)
			# for every other species' missing fragments, see if this fragment is in that species; assign presence if true
			for j in range(0, len(missing_target_frag_list)): 
				if i != j: # exclude self
					other_spec_frags = missing_target_frag_list[j]
					if frag in other_spec_frags: 
						presence.append(1)
					else: 
						presence.append(0) 
				elif i == j: 
						presence.append(1) # self: auto assign present
			oldest_present = presence.index(1) # finds leftmost 1 in array, ie oldest species in addition to the current species in which fragment is present
			if i != oldest_present: 
				continue # avoid double-counting fragments: this computes the branch assignment, so if you do it for all species for the same fragment, you will multi-count. 
					# arbitrarily do assignment when encounter it in the oldest species. 
			suspect_idx = [j for j in range(len(presence)) if j != i and presence[j] == 0 and coverage_rulings[j] == 0] # low-cov species marked not-shared (not-shared is the default when region can't be found)
			has_sharer = any(presence[j] == 1 for j in range(len(presence)) if j != i) # is there a confident sharer besides self
			if any(j < i for j in suspect_idx) or (not has_sharer and suspect_idx): # only exclude when a suspect could be the real oldest, or could flip single -> shared with no confident sharer
				continue
			spec_count = presence.count(1)
			oldest_other = next((idx for idx, val in enumerate(presence) if val == 1 and idx != i), None) # find the next oldest species
			youngest = len(presence) - 1 - presence[::-1].index(1) if 1 in presence else None # find the younget species
			if i == 0: # if oldest is helvetica
				if spec_count == 1: # only in helvetica
					branch = 0 # helvetica-only branch
				else: # if present in anyone else, assume in common ancestor
					branch = 8 # root
			elif i == 1: # if oldest is azteca
				if spec_count == 1: # only in azteca
					branch = 2 # azteca specific branch
				elif oldest_other == 0: # if in helvetica too, assume in common ancestor of both, which is 
					branch = 8 # root branch
				elif oldest_other == 2 or oldest_other == 3: # if in athabasca or algonquin too, assume in branch pre-their divergence
					branch = 1 # leading to dazt/dalg/daff/dath
			elif i == 2: # if oldest is algonquin
				if spec_count == 1:  # if only in algonquin
					branch = 6 # alg-specific branch
				elif oldest_other == 3: # if also in athabasca
					branch = 4 # in branch leading to dalg/dath
			elif i == 3: # if oldest is athabasca
				if spec_count == 1:  # if only in athabasca
					branch = 5 # athabasca-specific branch
			# add fragment to the selected branch in the dictionary
			region_branch_dict[region][branch].append(frag) 

branchnums = [0]*9 # initialize array that counts total number of fragments on each branch, from all regions
# 9 total branches including ancestral
# initialize array that has all fragments for each branch, from all regions
frags_by_branch = [[] for _ in range(numbranches+1)] 
# populate by summing over regions
for region in region_branch_dict: 
	branchfrags = region_branch_dict[region] 
	for i in range(0, len(branchfrags)): 
		branch = branchfrags[i]
		frags_by_branch[i].append(branch) # add list of fragments from region to master all-region branch frag list
		numfrags = len(branch) 
		branchnums[i] += numfrags

# total branch length in tree
totalbranchlen = sum(branchlens)

# test for overall heterogeneity across phylogeny per locus, but doesn't detect per branch; more power, less resolution
# (computed in memory; file output commented out)
#outfile2 = open('phylo_TL_analysis_fraggain_perlocus', 'w')
#outfile2.write('#Region \t Chi-square \t p-value \t Best_case_p-value \n')

omni_results = {} # region -> (omni_chi2, omni_p, best_p) or None if no flux
for region in region_branch_dict: 
	branchfrags = region_branch_dict[region] # fragments per branch in this TL
	totalvarfrags = 0 # track how many variable fragments there are at the locus
	for i in range(0, len(branchfrags)-1): # branch 8 is pre-clade (ancestral); no length, so can't do computation 
		branch = branchfrags[i]
		numfrags = len(branch) 
		totalvarfrags += numfrags
	# omnibus exact multinomial goodness-of-fit test for this locus
	if totalvarfrags > 0:
		observed = tuple(len(branchfrags[i]) for i in range(numbranches))
		branchprops = [branchlens[i]/totalbranchlen for i in range(numbranches)]
		omni_p, omni_chi2 = exact_multinomial_pvalue(observed, branchprops, totalvarfrags)
		shortest_branch = int(np.argmin(branchprops))
		best_obs = tuple(totalvarfrags if i == shortest_branch else 0 for i in range(numbranches))
		best_p, _ = exact_multinomial_pvalue(best_obs, branchprops, totalvarfrags)
		omni_results[region] = (omni_chi2, omni_p, best_p)
		#outfile2.write(region + '\t' + str(omni_chi2) + '\t' + str(omni_p) + '\t' + str(best_p) +  '\n')
	else: 
		omni_results[region] = None
		#outfile2.write(region + '\t' + '0' + '\t' + 'n/a' + '\t' + '1' + '\n') # best p-value is 1 so don't need to add another case when filter out high values here
#outfile2.close()

## per-branch, per-locus binomial tests
outfile = open('Phylogenetic_analysis_fragment_gain_per_branch_per_locus', 'w')
outfile.write('#Region \t Branch \t Expected_frags \t Observed_frags \t p-value \t Best_case_p-value \n')

for region in region_branch_dict: 
	branchfrags = region_branch_dict[region]
	totalvarfrags = sum(len(branchfrags[i]) for i in range(numbranches))
	for i in range(0, len(branchfrags)-1): # branch 8 is pre-clade (ancestral); no length
		if totalvarfrags == 0: 
			outfile.write(region + '\t' + 'Branch_' + str(i) + '\t' + str(0) + '\t' + str(0) + '\t' + 'n/a' + '\t' + '1' + '\n')
			continue
		numfrags = len(branchfrags[i])
		proplen = branchlens[i]/totalbranchlen
		expnumfrags = proplen*totalvarfrags
		prob = branchlens[i]/totalbranchlen
		p = binomtest(numfrags, totalvarfrags, prob, alternative='two-sided').pvalue
		bestp = binomtest(totalvarfrags, totalvarfrags, prob, alternative='two-sided').pvalue
		outfile.write(region + '\t' + 'Branch_' + str(i) + '\t' + str(expnumfrags) + '\t' + str(numfrags) + '\t' + str(p) + '\t' + str(bestp) + '\n')
outfile.close()

## now do an analysis summing over all fragments in all loci
# to see if rate of total fragment gain is different aross branches

# initialize list of total fragment gains for each branch
branch_frag_counts = [0]*(numbranches) #branch 8 is pre-clade (ancestral); no length
# sum across loci, excepting branch 8
for region in region_branch_dict: 
	branchfrags = region_branch_dict[region] 
	for i in range(0, len(branchfrags)-1): # branch 8 is pre-clade (ancestral); no length
		branch = branchfrags[i]
		numfrags = len(branch) 
		branch_frag_counts[i] += numfrags 

# do binomial calculation on each branch and print
outfile = open('Phylogenetic_analysis_fragment_gain_aggregated_per_branch', 'w')
outfile.write('#Branch \t Expected_frags \t Observed_frags \t p-value \t Best_case_p-value \n')

totalvarfrags = np.sum(branch_frag_counts) 
for i in range(0, numbranches): 
	numfrags = branch_frag_counts[i]
	prob = branchlens[i]/totalbranchlen
	expnumfrags = prob*totalvarfrags
	p = binomtest(numfrags, totalvarfrags, prob, alternative='two-sided').pvalue
	bestp = binomtest(totalvarfrags, totalvarfrags, prob, alternative='two-sided').pvalue
	outfile.write('Branch_' + str(i) + '\t' + str(expnumfrags) + '\t' + str(numfrags) + '\t' + str(p) + '\t' + str(bestp) + '\n')
outfile.close()

# now go back to raw p value output files for each of the 3 files and do FDR correction for all of the entities (branches, loci) that have sufficient power

alpha = 0.05 # use same alpha for FDR and filtering

## first, aggregated fragments over loci per branch
infile = np.genfromtxt('Phylogenetic_analysis_fragment_gain_aggregated_per_branch', dtype=str, delimiter='\t', skip_header=True)
outfile = open('Phylogenetic_analysis_fragment_gain_aggregated_per_branch_FDR', 'w')
outfile.write('#Branch \t Expected_frags \t Observed_frags \t p-value(Corrected) \t Significant? \n')


sigcount = 0
nonsigcount = 0
nopowercount = 0
# filter rows based on having power to count total number of effective tests for p-value correction
filtidx = []
for i in range(len(infile)): 
	row = infile[i]
	bestp = float(row[4]) 
	if bestp < alpha: 
		filtidx.append(i) 

# get filtered p values for FDR
filtpvals = [float(infile[i][3]) for i in filtidx]
# perform fdr
_, qvals, _, _ = multipletests(filtpvals, alpha=alpha, method='fdr_bh')
qval_lookup = dict(zip(filtidx, qvals))

for i in range(len(infile)): 
	row = infile[i]
	if i in filtidx:
		qval = qval_lookup[i]
		row[3] = str(qval) # replace p value with FDR corrected p value
		if qval < alpha: 
			outfile.write('\t'.join(row[:-1]) + '\t' + 'TRUE' + '\n') # write row with p-value replaced and the best case p value column replaced with indicator of significant or of lack of power
			sigcount += 1
		else: 
			outfile.write('\t'.join(row[:-1]) + '\t' + 'FALSE' + '\n')
			nonsigcount += 1
	else: 
		outfile.write('\t'.join(row[:-1]) + '\t' + 'N/A' + '\n')
		nopowercount += 1
outfile.close()

print('Phylogenetic_analysis_fragment_gain_aggregated_per_branch')
print('no power:', nopowercount)
print('significant:', sigcount) 
print('not significant', nonsigcount)

## next, omnibus chi square test
# (FDR computation and file output commented out)
#infile = np.genfromtxt('phylo_TL_analysis_fraggain_perlocus', dtype=str, delimiter='\t', skip_header=True)
#outfile = open('phylo_TL_analysis_fraggain_perlocus_FDR', 'w')
#outfile.write('#Region \t Chi-square \t p-value(Corrected) \t Significant? \n')
#
#sigcount = 0
#nonsigcount = 0
#nopowercount = 0
## filter rows based on having power to count total number of effective tests for p-value correction
#filtidx = []
#for i in range(len(infile)): 
#        row = infile[i]
#        bestp = float(row[3])  # fourth column here
#        if bestp < alpha: 
#                filtidx.append(i)
# 
## get filtered p values for FDR
#filtpvals = [float(infile[i][2]) for i in filtidx]
## perform fdr
#_, qvals, _, _ = multipletests(filtpvals, alpha=alpha, method='fdr_tsbh')
#qval_lookup = dict(zip(filtidx, qvals))
#
#for i in range(len(infile)): 
#	row = infile[i]
#	if i in filtidx:
#		qval = qval_lookup[i]
#		row[2] = str(qval)
#		if qval < alpha: 
#			outfile.write('\t'.join(row[:-1]) + '\t' + 'TRUE' + '\n') # write row with p-value replaced and the best case p value column replaced with indicator of significant or of lack of power
#			sigcount += 1
#		else: 
#			outfile.write('\t'.join(row[:-1]) + '\t' + 'FALSE' + '\n')
#			nonsigcount += 1
#	else: 
#		outfile.write('\t'.join(row[:-1]) + '\t' + 'N/A' + '\n')
#		nopowercount += 1
#outfile.close()
#
#print('phylo_TL_analysis_fraggain_perlocus')
#print('no power:', nopowercount)
#print('significant:', sigcount) 
#print('not significant', nonsigcount)


## enrichment of low p-values
# for each threshold T: count how many tests have raw p < T, compare to expected count under H0 (exact, computed from the discrete null distribution)

thresholds = [0.05]
branchprops = [branchlens[i]/totalbranchlen for i in range(numbranches)]

# (b) omnibus per-locus exact multinomial tests
# (computation kept; print commented out)
omni = []
for region in region_branch_dict: 
	branchfrags = region_branch_dict[region]
	n = sum(len(branchfrags[i]) for i in range(numbranches))
	if n == 0: 
		continue
	observed_counts = tuple(len(branchfrags[i]) for i in range(numbranches))
	rawp, _ = exact_multinomial_pvalue(observed_counts, branchprops, n)
	shortest = int(np.argmin(branchprops))
	best_obs = tuple(n if i == shortest else 0 for i in range(numbranches))
	bestp, _ = exact_multinomial_pvalue(best_obs, branchprops, n)
	omni.append((rawp, n, bestp))

# (c) locus-level: how many loci have at least one branch with p < T?
# expected count computed under exact joint multinomial null, NOT independence-union of marginal binomials
print('\nenrichment of loci with at least one significant branch (locus-level test):')
locus_branches = {}
for region in region_branch_dict: 
	branchfrags = region_branch_dict[region]
	n = sum(len(branchfrags[i]) for i in range(numbranches))
	if n == 0: 
		continue
	locus_branches[region] = []
	for i in range(numbranches): 
		prob = branchlens[i]/totalbranchlen
		k = len(branchfrags[i])
		rawp = binomtest(k, n, prob, alternative='two-sided').pvalue
		bestp = binomtest(n, n, prob, alternative='two-sided').pvalue
		locus_branches[region].append((rawp, n, prob, bestp))

for T in thresholds: 
	observed = 0
	expected = 0.0
	m_capable = 0
	for region, branches in locus_branches.items(): 
		n = branches[0][1] # same n across branches within a locus
		probs_locus = [prob for (rp, n_, prob, bp) in branches]
		bps = [bp for (rp, n_, prob, bp) in branches]
		rps = [rp for (rp, n_, prob, bp) in branches]
		# capable if at least one branch could possibly yield p<T at this n
		if not any(bp < T for bp in bps): 
			continue
		m_capable += 1
		if any(rp < T for rp in rps): 
			observed += 1
		# precompute the set of k values that give p<T for each branch (branch-specific because prob differs)
		sig_k = []
		for i in range(numbranches): 
			s = set()
			for k in range(n+1): 
				if binomtest(k, n, probs_locus[i], alternative='two-sided').pvalue < T: 
					s.add(k)
			sig_k.append(s)
		# enumerate multinomial compositions and accumulate joint P(any branch p<T)
		logps = []
		indicators = []
		for comp in enumerate_compositions(n, numbranches): 
			logps.append(log_multinomial_pmf(comp, probs_locus, n))
			ind = 1 if any(comp[i] in sig_k[i] for i in range(numbranches)) else 0
			indicators.append(ind)
		logps_arr = np.array(logps)
		ind_arr = np.array(indicators)
		mx = logps_arr.max()
		w = np.exp(logps_arr - mx)
		p_any = float((w * ind_arr).sum() / w.sum())
		expected += p_any
	excess = observed - expected
	# two-tailed Poisson test: take one-tailed p in the direction of deviation, double it (capped at 1)
	if observed >= expected: 
		one_tail = 1 - poisson.cdf(observed - 1, expected) if observed > 0 else 1.0
	else: 
		one_tail = poisson.cdf(observed, expected)
	enrich_p = min(2 * one_tail, 1.0)
	print('  threshold p<' + str(T) + ': m_capable=' + str(m_capable) + ', observed=' + str(observed) + ', expected=' + str(round(expected, 2)) + ', excess=' + str(round(excess, 2)) + ', two-tailed p=' + str(enrich_p))

# per-branch counts: how many loci had each branch's raw p below T?
print('\nper-branch counts of loci with p < T:')
for T in thresholds: 
	counts = [0]*numbranches
	for region, branches in locus_branches.items(): 
		for i, (rp, n, prob, bp) in enumerate(branches): 
			if rp < T: 
				counts[i] += 1
	print('  threshold p<' + str(T) + ':', counts)
