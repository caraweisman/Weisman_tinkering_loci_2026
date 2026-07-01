mport numpy as np
from statsmodels.stats.multitest import multipletests
import sys
import scipy
from scipy.stats import binomtest, binom, poisson
import re 

## ***** counts each fragment, even if in same gene (must, because these are duplicates) 
# assigns gain of fragment(s) from each unique gene to a branch in the phylogeny
# does a statistical test to see if a given locus has a nonuniform rate of fragment gain across the phylogeny
# updates 5/25: - added chi square test for overall test for heterogeniety; above test is per-branch, this is not. higher power, lower resolution
# - changed from bonferroni correction here to second sweep to FDR correct given the best case p-values as the number of tests
# updated 5/30: added new filter for removing ambiguous fragments. test whether their ambiguity matters for the result, and if not keep; if so, exclude.
# updated: locus-level "at least one significant branch" enrichment test now uses joint multinomial null (MC) instead of independence-assumption union formula


species = ['Dhel', 'Dazt', 'Dalg', 'Dath']
speciestree = [0,1,2,3]

regions = np.genfromtxt('TL_sfetch', dtype=str, delimiter='\t', usecols=(0), comments='#')

resultfiles = ['duplicate_flux_' + s for s in species]


def long_to_short_name(name): # flux output files include framgent locations; remove this in conservative assumption that they all stem from the same gain event
        desc = '-'.join(re.split('-', name)[3:])
        isofree = re.split(',', desc)[0]
        return isofree


from scipy.stats import multinomial as multinomial_dist

def chisq_stat(counts, expected): 
	s = 0.0
	for c, e in zip(counts, expected): 
		if e > 0: 
			s += (c - e)**2 / e
	return s

def exact_multinomial_pvalue(observed, probs, n, n_resamples=99999): 
	# Monte Carlo estimate of multinomial GOF p-value (drop-in replacement for exact enumeration)
	expected = [n * p for p in probs]
	obs_stat = chisq_stat(observed, expected)
	draws = multinomial_dist.rvs(int(n), probs, size=n_resamples)
	# vectorized chi-square statistic across resamples
	exp_arr = np.array(expected)
	null_stats = ((draws - exp_arr)**2 / np.where(exp_arr > 0, exp_arr, 1)).sum(axis=1)
	# p-value = (number of null stats >= observed + 1) / (n_resamples + 1) (standard MC correction)
	pval = (1 + (null_stats >= obs_stat - 1e-12).sum()) / (n_resamples + 1)
	return float(pval), obs_stat


# Sankoff DP algorithm to compute losses/gains per branch by parsimony
# In the case of a tie, prefers the ancestral branch
# Two versions to account for the difference in meaning of the flux values
# optional agnostic flag to compare results with and without agnostic leaves due to low coverage and exclude cases where they conflict

def parsimony_events_focal(counts, agnostic=()):
	"""
	counts: [Dhel, Dazt, Dalg, Dath] - missing counts per species relative to affinis.
	agnostic: set of species names whose leaves should be treated as unknown.
	Returns: list of 9 ints, events per branch.
	"""
	aff_count = max(counts)
	leaf_counts = {
		'Dhel': aff_count - counts[0],
		'Dazt': aff_count - counts[1],
		'Dalg': aff_count - counts[2],
		'Dath': aff_count - counts[3],
		'Daff': aff_count,
	}
	return _run_parsimony(leaf_counts, agnostic)
 
 
def parsimony_events_target(counts, agnostic=()):
	"""
	counts: [Dhel, Dazt, Dalg, Dath] - target-side fragment counts per species (not in affinis).
	agnostic: set of species names whose leaves should be treated as unknown.
	Returns: list of 9 ints, events per branch.
	"""
	leaf_counts = {
		'Dhel': counts[0],
		'Dazt': counts[1],
		'Dalg': counts[2],
		'Dath': counts[3],
		'Daff': 0,
	}
	return _run_parsimony(leaf_counts, agnostic)
 
 
def _run_parsimony(leaf_counts, agnostic=()):
	"""Shared Sankoff DP + traceback. Takes a dict of leaf names -> counts."""
	INF = float('inf')
	K = max(leaf_counts.values()) + 2
 
	def leaf_S(name):
		if name in agnostic:
			return [0] * K                                  # uninformative: any count consistent
		c = leaf_counts[name]
		return [0 if k == c else INF for k in range(K)]
 
	def internal_S(child_S_list):
		S = []
		for k in range(K):
			total = 0
			for cs in child_S_list:
				total += min(cs[j] + abs(k - j) for j in range(K))
			S.append(total)
		return S
 
	S_Dhel = leaf_S('Dhel')
	S_Dazt = leaf_S('Dazt')
	S_Daff = leaf_S('Daff')
	S_Dalg = leaf_S('Dalg')
	S_Dath = leaf_S('Dath')
	S_n3 = internal_S([S_Dalg, S_Dath])
	S_n2 = internal_S([S_Daff, S_n3])
	S_n1 = internal_S([S_Dazt, S_n2])
	S_root = internal_S([S_Dhel, S_n1])
 
	root_count = min(range(K), key=lambda k: (S_root[k] + k, -k))
 
	def best_child_count(S_child, parent_count):
		return min(range(K), key=lambda j: S_child[j] + abs(parent_count - j))
 
	n1_count = best_child_count(S_n1, root_count)
	n2_count = best_child_count(S_n2, n1_count)
	n3_count = best_child_count(S_n3, n2_count)
	# resolve leaf counts too: a normal leaf returns its observed count,
	# an agnostic leaf snaps to its parent -> zero events on that branch
	dhel = best_child_count(S_Dhel, root_count)
	dazt = best_child_count(S_Dazt, n1_count)
	daff = best_child_count(S_Daff, n2_count)
	dalg = best_child_count(S_Dalg, n3_count)
	dath = best_child_count(S_Dath, n3_count)
 
	events = [0] * 9
	events[0] = abs(root_count - dhel)
	events[1] = abs(root_count - n1_count)
	events[2] = abs(n1_count - dazt)
	events[3] = abs(n1_count - n2_count)
	events[4] = abs(n2_count - n3_count)
	events[5] = abs(n3_count - dalg)
	events[6] = abs(n3_count - dath)
	events[7] = abs(n2_count - daff)
	events[8] = root_count
	return events
 
 

# initialize dictionary that will hold compiled results from all species
region_flux_dict = {}
for region in regions: 
	region_flux_dict[region] = [] # for every region, an entry in the dict
	for s in species: 
		region_flux_dict[region].append([[],[]]) 
# for every region, the dict entry gets one array per species, consisting of two subarrays: things present in focal but missing
# in this species, and present in this species and missin from focal


# populate dictionary with data from individual species' flux filesspeciesidx = 0  
speciesidx = 0 
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
					region_flux_dict[region][speciesidx][0].append(long_to_short_name(hit)) 
			else: 
				for hit in hits: 
					region_flux_dict[region][speciesidx][1].append(long_to_short_name(hit)) 
	speciesidx += 1

# assign branches and lengths
numbranches = 8 # 8 within the phylogeny with lengths; ancestral branch (8) is of indeterminate length
#branchlens = [3.9,1.8,2.1,0.6,0.2,1.3,1.3,1.5] # from BK HOG UCLD tree, PLOS Biology 2024 (is this right?)
branchlens = [0.063, 0.041, 0.063, 0.013, 0.005, 0.028, 0.035, 0.046] # from Kim PLOS biology 2024 actually - 

# initialize dictionary that will hold fragment assignments per branch per region, to be populated below in loop
region_branch_dict = {}

# load coverage information for each region
# only applies to affinis; if cov < threshold in a species, don't consider
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
	# get coverage for each region/species
	coverage_rulings = coverage_dict[region]
	# initialize dictionary that will hold the fragments assigned to each branch for this region
	region_branch_dict[region] = [[] for _ in range(numbranches + 1)]  # +1 for the ancestral branch

	# focal-to-target
	all_focal_frag_names = {f for i in region_flux_dict[region] for f in i[0]}
	missing_aff_frag_dict = {}
	for frag in all_focal_frag_names:
		counts = [region_flux_dict[region][s][0].count(frag) for s in range(len(species))]
		# low-cov species that look like they match affinis (counts == 0) are the suspect/default calls
		suspect = {species[s] for s in range(len(species)) if coverage_rulings[s] == 0 and counts[s] == 0}
		events_face = parsimony_events_focal(counts)            # trust the zeros
		events_agn = parsimony_events_focal(counts, suspect)    # treat the suspect zeros as unknown
		if events_face != events_agn:
			continue                                        # answer hinges on untrusted data -> exclude
		missing_aff_frag_dict[frag] = counts
		branch_events = events_face                             # == events_agn here
		for i in range(0, len(region_branch_dict[region])):
			for _ in range(branch_events[i]):
				region_branch_dict[region][i].append(frag)

	# target-to-focal
	all_target_frag_names = {f for i in region_flux_dict[region] for f in i[1]}
	missing_target_frag_dict = {}
	for frag in all_target_frag_names:
		counts = [region_flux_dict[region][s][1].count(frag) for s in range(len(species))]
		suspect = {species[s] for s in range(len(species)) if coverage_rulings[s] == 0 and counts[s] == 0}
		events_face = parsimony_events_target(counts)
		events_agn = parsimony_events_target(counts, suspect)
		if events_face != events_agn:
			continue
		missing_target_frag_dict[frag] = counts
		branch_events = events_face
		for i in range(0, len(region_branch_dict[region])):
			for _ in range(branch_events[i]):
				region_branch_dict[region][i].append(frag)

branchnums = [0]*9 # initialize array that counts total number of fragments on each branch, from all regions
# 9 total branches including ancestral
# initialize array that has all fragments for each branch, from all regions
frags_by_branch = [[] for _ in range(numbranches+1)]  
# populate by summing over regions
for region in region_branch_dict: 
	branchfrags = region_branch_dict[region] 
	for i in range(0, len(branchfrags)): 
		branch = branchfrags[i]
		frags_by_branch[i].append(branch)
		numfrags = len(branch) 
		branchnums[i] += numfrags

# total branch length in tree
totalbranchlen = sum(branchlens)

# test for overall heterogeneity across phylogeny per locus, but doesn't detect per branch; more power, less resolution
# (computed in memory; file output commented out)
#outfile2 = open('phylo_TL_analysis_fragdup_perlocus', 'w')
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
		#outfile2.write(region + '\t' + '0' + '\t' + 'n/a' + '\t' + '1' + '\n')
#outfile2.close()

## per-branch, per-locus binomial tests
outfile = open('Phylogenetic_analysis_fragment_duplication_per_branch_per_locus', 'w')
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
outfile = open('Phylogenetic_analysis_fragment_duplication_aggregated_per_branch', 'w')
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

# now go back to raw p value output files and do FDR correction for all of the entities (branches, loci) that have sufficient power

alpha = 0.05 # use same alpha for FDR and filtering

## first, aggregated fragments over loci per branch
infile = np.genfromtxt('Phylogenetic_analysis_fragment_duplication_aggregated_per_branch', dtype=str, delimiter='\t', skip_header=True)
outfile = open('Phylogenetic_analysis_fragment_duplication_aggregated_per_branch_FDR', 'w')
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

print('Phylogenetic_analysis_fragment_duplication_aggregated_per_branch')
print('no power:', nopowercount)
print('significant:', sigcount) 
print('not significant', nonsigcount)

## next, omnibus chi square test
# (FDR computation and file output commented out)
#infile = np.genfromtxt('phylo_TL_analysis_fragdup_perlocus', dtype=str, delimiter='\t', skip_header=True)
#outfile = open('phylo_TL_analysis_fragdup_perlocus_FDR', 'w')
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
#print('phylo_TL_analysis_fragdup_perlocus')
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

#print('\nenrichment of low p-values: omnibus per-locus')
#for T in thresholds: 
#	capable = [(rp, n, bp) for (rp, n, bp) in omni if bp < T]
#	observed = sum(1 for (rp, n, bp) in capable if rp < T)
#	# expected: for each capable locus, P(p<T|H0) estimated by Monte Carlo
#	# draw multinomial samples, compute each draw's chi-square stat, compute each draw's own p-value (rank among null), count fraction with p < T
#	expected = 0.0
#	for (rp, n, bp) in capable: 
#		exp_counts = np.array([n * p for p in branchprops])
#		draws = multinomial_dist.rvs(int(n), branchprops, size=99999)
#		null_stats = ((draws - exp_counts)**2 / np.where(exp_counts > 0, exp_counts, 1)).sum(axis=1)
#		# each draw's own p-value = fraction of null_stats >= its stat (MC estimate). vectorize via rank.
#		# count outcomes whose p-value < T: equivalent to counting outcomes in the top T fraction of null_stats by rank
#		# (with the MC correction (1+k)/(N+1), p < T iff rank from top < T*(N+1) - 1, approx top T fraction)
#		num_with_p_below_T = int(T * (len(null_stats) + 1))
#		# get threshold stat at the top T fraction
#		if num_with_p_below_T <= 0: 
#			frac_below = 0.0
#		else: 
#			threshold_stat = np.partition(null_stats, -num_with_p_below_T)[-num_with_p_below_T]
#			frac_below = (null_stats >= threshold_stat).sum() / len(null_stats)
#		expected += frac_below
#	excess = observed - expected
#	enrich_p = 1 - poisson.cdf(observed - 1, expected) if observed > 0 else 1.0
#	print('  threshold p<' + str(T) + ': m_capable=' + str(len(capable)) + ', observed=' + str(observed) + ', expected=' + str(round(expected, 2)) + ', excess=' + str(round(excess, 2)) + ', enrichment p=' + str(enrich_p))

# (c) locus-level: how many loci have at least one branch with p < T?
# expected count computed under joint multinomial null via MC, NOT independence-union of marginal binomials
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
	n_mc = 99999
	for region, branches in locus_branches.items(): 
		n = branches[0][1]
		probs_locus = [prob for (rp, n_, prob, bp) in branches]
		bps = [bp for (rp, n_, prob, bp) in branches]
		rps = [rp for (rp, n_, prob, bp) in branches]
		if not any(bp < T for bp in bps): 
			continue
		m_capable += 1
		if any(rp < T for rp in rps): 
			observed += 1
		# MC estimate of P(any branch p<T) under joint multinomial null
		# precompute the set of k values that give p<T for each branch (branch-specific because prob differs)
		sig_k = []
		for i in range(numbranches): 
			s = set()
			for k in range(n+1): 
				if binomtest(k, n, probs_locus[i], alternative='two-sided').pvalue < T: 
					s.add(k)
			sig_k.append(s)
		draws = multinomial_dist.rvs(int(n), probs_locus, size=n_mc)
		# per draw: 1 if any branch's count is in that branch's sig_k set
		any_sig = np.zeros(n_mc, dtype=bool)
		for i in range(numbranches): 
			if sig_k[i]: 
				any_sig |= np.isin(draws[:, i], list(sig_k[i]))
		p_any = (1 + any_sig.sum()) / (n_mc + 1) # standard MC correction
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
