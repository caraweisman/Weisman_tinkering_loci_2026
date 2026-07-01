import numpy as np
import sys 
import re 
from scipy.stats import binomtest
from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact
import random 
from statsmodels.stats.multitest import multipletests

print('starting', flush=True)

# usage: 7recomb_site_enrichment.py (TE gff - one from annotation) (gff type - one of DFAM, HMMER, RNAPROT, RNAREAD) (contig lengths) (flank) (main gff) (TL gff)

# flank TE density analysis only; no recombination tracts
# four groups of features, each feature = two flank regions of 'flank' nt:
#   a) nonparental conserved exons: conserved exons that do NOT match any fragment
#   b) parental conserved exons: conserved exons that DO match a fragment (same accession, >0 aa overlap)
#   c) fragments
#   d) intergenic controls: pretend-fragments placed in clean intergenic space; body length copied from each real fragment, flanked the same way; only the flanks are measured
# parental matching uses the gffs directly (accession + aa coordinate overlap), NOT a recomb_hits file
# overall (any TE) enrichment: mann-whitney on the per-feature percent-covered distributions; two comparisons (nonparental vs parental, fragments vs controls); BH over the 2 tests
# per-family enrichment: same two comparisons; mann-whitney per family on the per-feature percent-covered values, and fisher exact per family on presence/absence per feature (covered > 0); BH over 2 * n_families, separately for the mwu set and the fisher set

TE_gff = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t') 
gff_type = sys.argv[2]
flank = int(sys.argv[4])
main_gff = np.genfromtxt(sys.argv[5], dtype=str, delimiter='\t', comments=None)
TL_gff = np.genfromtxt(sys.argv[6], dtype=str, delimiter='\t', comments=None)

# contig lengths, for sampling control regions
contiglenfile = np.genfromtxt(sys.argv[3], dtype=str)
contig_lengths = {}
for line in contiglenfile: 
	contig = line[0]
	length = int(line[1])
	contig_lengths[contig] = length

def extract_id_dfam(desc): 
	teid = re.split('-', re.split('_', desc)[0])[-2:]
	return '-'.join(teid)

def extract_id_hmmer(desc): 
	teid = re.split('_JA', desc)[0][3:]
	return teid

# conserved exons: flank regions for coverage, plus parsed accession + aa-coords for parental matching
# consname = accession + '_' + aa-coords, eg NP_001260771.1_278-382; one entry per exon
main_gff_dict = {} # key = consname, value = [us, ds] flank regions
cons_exons = [] # one row per exon: [accession, aa_start, aa_stop, consname]
for row in main_gff: 
	desc = row[8]
	acc = re.split('-', desc)[2]
	coords = re.split('/', desc)[0][4:]
	consname = acc + '_' + coords
	contig = row[0]
	start = int(row[3])
	stop = int(row[4]) 
	us = [contig, start-flank, start]
	ds = [contig, stop, stop+flank]
	main_gff_dict[consname] = []
	main_gff_dict[consname].append(us)
	main_gff_dict[consname].append(ds)
	aa_start = int(re.split('-', coords)[0])
	aa_stop = int(re.split('-', coords)[1])
	cons_exons.append([acc, aa_start, aa_stop, consname])

TE_pos_dict = {} # dictionary for all positions of TEs
for row in TE_gff: 
	contig = row[0]
	start = int(row[3])
	stop = int(row[4]) 
	if gff_type == 'DFAM': 
		teid = extract_id_dfam(row[8]) 
	elif gff_type == 'HMMER':
		teid = extract_id_hmmer(row[8])
	elif gff_type == 'RNAPROT': 
		teid = row[8]
	if contig not in TE_pos_dict: 
		TE_pos_dict[contig] = []
	TE_pos_dict[contig].append([start, stop, teid]) 

# fragments: parse accession + aa-coords for parental matching, and record flank regions + body length
frag_exons = [] # one row per fragment: [accession, aa_start, aa_stop]
frag_regions = [] # one row per fragment: [us, ds] flank regions
frag_body_lengths = [] # one body length per fragment, for building controls
for row in TL_gff: 
	desc = row[8]
	acc = re.split('-', desc)[2]
	coords = re.split('/', desc)[0][4:]
	aa_start = int(re.split('-', coords)[0])
	aa_stop = int(re.split('-', coords)[1])
	contig = row[0]
	start = int(row[3])
	stop = int(row[4])
	frag_exons.append([acc, aa_start, aa_stop, contig, start, stop])
	us = [contig, start-flank, start]
	ds = [contig, stop, stop+flank]
	frag_regions.append([us, ds])
	frag_body_lengths.append(stop - start)

# parental exons: a conserved exon is parental if some fragment has the same accession and >0 aa overlap
# flat double loop over fragments x conserved exons; mark the consname
parental = {} # key = consname of a parental conserved exon, value = 1
n_unmatched = 0 # fragments with no matching conserved exon
for frag in frag_exons: 
	frag_acc = frag[0]
	frag_aa_start = frag[1]
	frag_aa_stop = frag[2]
	frag_contig = frag[3]
	frag_start = frag[4]
	frag_stop = frag[5]
	for exon in cons_exons: 
		exon_acc = exon[0]
		exon_aa_start = exon[1]
		exon_aa_stop = exon[2]
		consname = exon[3]
		if exon_acc != frag_acc: 
			continue
		# >0 aa overlap: the two aa ranges share at least one position
		if max(frag_aa_start, exon_aa_start) <= min(frag_aa_stop, exon_aa_stop): 
			parental[consname] = 1

# split the conserved exons into nonparental and parental (nonoverlapping sets)
nonparental_feats = [] # each feature is a list of flank regions [us, ds]
parental_feats = []
for consname in main_gff_dict: 
	regions = main_gff_dict[consname]
	if consname in parental: 
		parental_feats.append(regions)
	else: 
		nonparental_feats.append(regions)

# control sampling: forbidden intervals are the bodies and both flanks of every main exon and every fragment
# a candidate control must not overlap any of these on its contig
forbidden = {} # contig -> list of [lo, hi]
for row in main_gff: 
	contig = row[0]
	start = int(row[3])
	stop = int(row[4])
	if contig not in forbidden: 
		forbidden[contig] = []
	forbidden[contig].append([start, stop]) # body
	forbidden[contig].append([start-flank, start]) # us flank
	forbidden[contig].append([stop, stop+flank]) # ds flank
for row in TL_gff: 
	contig = row[0]
	start = int(row[3])
	stop = int(row[4])
	if contig not in forbidden: 
		forbidden[contig] = []
	forbidden[contig].append([start, stop]) # body
	forbidden[contig].append([start-flank, start]) # us flank
	forbidden[contig].append([stop, stop+flank]) # ds flank

# for sampling a contig in proportion to its length
contig_names = []
contig_weights = []
for contig in contig_lengths: 
	contig_names.append(contig)
	contig_weights.append(contig_lengths[contig])

# intergenic controls: one per fragment, body length copied from that fragment
# place a pretend-fragment of length L in clean intergenic space and take its two flanks; only the flanks are measured, so the structure matches a real fragment exactly (two 'flank' nt regions separated by a body of length L)
max_attempts = 100000
control_feats = []
for L in frag_body_lengths: 
	placed = False
	attempts = 0
	while placed == False and attempts < max_attempts: 
		attempts += 1
		contig = random.choices(contig_names, weights=contig_weights)[0]
		clen = contig_lengths[contig]
		# need room for the upstream flank, the body, and the downstream flank
		lowest = flank
		highest = clen - L - flank
		if highest < lowest: 
			continue # contig too short for this length; try another
		s = random.randint(lowest, highest)
		span_lo = s - flank
		span_hi = s + L + flank
		# reject if the whole span overlaps anything forbidden on this contig
		bad = False
		contig_forbidden = forbidden.get(contig, [])
		for interval in contig_forbidden: 
			ilo = interval[0]
			ihi = interval[1]
			if max(span_lo, ilo) <= min(span_hi, ihi): 
				bad = True
				break
		if bad == True: 
			continue
		us = [contig, s-flank, s]
		ds = [contig, s+L, s+L+flank]
		control_feats.append([us, ds])
		placed = True
	if placed == False: 
		print('warning: could not place a control region for body length', L, 'after', max_attempts, 'attempts', flush=True)

# TE density in flanking regions (nt annotated as TE / 2*FLANK nt) computed per feature, for each of the four groups
# for each feature we record both the overall covered fraction (any TE) and the per-family covered fraction (one value per family)
groups = [] # each entry: [name, list of features]; each feature is a list of flank regions
groups.append(['nonparental_exons', nonparental_feats])
groups.append(['parental_exons', parental_feats])
groups.append(['fragments', frag_regions])
groups.append(['intergenic_controls', control_feats])

group_fracs = {} # name -> list of overall percent-covered values, one per feature
group_fam_frac = {} # name -> list of per-feature dicts (TEid -> percent-covered by that family)
for group in groups: 
	name = group[0]
	feats = group[1]
	fracs = []
	feat_fam_fracs = []
	for regions in feats: 
		covered = 0
		fam_covered = {} # key = TEid, value = total covered bases by that family over this feature's two flank regions
		for region in regions: 
			contig = region[0]
			regstart = region[1]
			regstop = region[2]
			marked = [False] * (regstop - regstart)
			fam_marked = {} # key = TEid, value = boolean array over this region
			contigTEs = TE_pos_dict.get(contig, [])
			for TE in contigTEs: 
				TEstart = TE[0]
				TEstop = TE[1]
				TEid = TE[2]
				ostart = max(TEstart, regstart)
				ostop = min(TEstop, regstop)
				if ostop >= ostart: 
					if TEid not in fam_marked: 
						fam_marked[TEid] = [False] * (regstop - regstart)
					for pos in range(ostart, ostop + 1): 
						idx = pos - regstart
						if idx >= 0 and idx < regstop - regstart: 
							marked[idx] = True
							fam_marked[TEid][idx] = True
			for k in range(0, regstop - regstart): 
				if marked[k] == True: 
					covered += 1
			for TEid in fam_marked: 
				famcov = 0
				for k in range(0, regstop - regstart): 
					if fam_marked[TEid][k] == True: 
						famcov += 1
				if TEid not in fam_covered: 
					fam_covered[TEid] = 0
				fam_covered[TEid] += famcov
		fracs.append(covered / float(2 * flank))
		feat_fam_frac = {}
		for TEid in fam_covered: 
			feat_fam_frac[TEid] = fam_covered[TEid] / float(2 * flank)
		feat_fam_fracs.append(feat_fam_frac)
	group_fracs[name] = fracs
	group_fam_frac[name] = feat_fam_fracs

nonparental_fracs = group_fracs['nonparental_exons']
parental_fracs = group_fracs['parental_exons']
fragment_fracs = group_fracs['fragments']
control_fracs = group_fracs['intergenic_controls']

# overall (any TE) enrichment: mann-whitney U tests (two-sided) for the two comparisons, BH over the 2 tests
u_np, p_np = mannwhitneyu(nonparental_fracs, parental_fracs, alternative='two-sided')
u_fc, p_fc = mannwhitneyu(fragment_fracs, control_fracs, alternative='two-sided')

density_pvals = [p_np, p_fc]
density_qvals = multipletests(density_pvals, method='fdr_bh')[1]

# direction of each overall comparison, from the U statistic (U is for the first group; U > n1*n2/2 means the first group has the larger values)
if u_np > (len(nonparental_fracs) * len(parental_fracs)) / 2.0: 
	direction_np = 'nonparental higher'
elif u_np < (len(nonparental_fracs) * len(parental_fracs)) / 2.0: 
	direction_np = 'parental higher'
else: 
	direction_np = 'equal'
if u_fc > (len(fragment_fracs) * len(control_fracs)) / 2.0: 
	direction_fc = 'fragments higher'
elif u_fc < (len(fragment_fracs) * len(control_fracs)) / 2.0: 
	direction_fc = 'controls higher'
else: 
	direction_fc = 'equal'

print('flank TE density - nonparental conserved exons: n', len(nonparental_fracs), 'mean', np.mean(nonparental_fracs), 'median', np.median(nonparental_fracs))
print('flank TE density - parental conserved exons: n', len(parental_fracs), 'mean', np.mean(parental_fracs), 'median', np.median(parental_fracs))
print('flank TE density - fragments: n', len(fragment_fracs), 'mean', np.mean(fragment_fracs), 'median', np.median(fragment_fracs))
print('flank TE density - intergenic controls: n', len(control_fracs), 'mean', np.mean(control_fracs), 'median', np.median(control_fracs))
print('MWU nonparental vs parental exons: p', p_np, 'q', density_qvals[0], 'direction', direction_np)
print('MWU fragments vs intergenic controls: p', p_fc, 'q', density_qvals[1], 'direction', direction_fc)

outfile2 = open('nonparent_exons_TE_density_' + gff_type, 'w')
outfile3 = open('parent_exons_TE_density_' + gff_type, 'w')
outfile4 = open('fragment_TE_density_' + gff_type, 'w')
outfile6 = open('intergenic_control_TE_density_' + gff_type, 'w')

for num in nonparental_fracs: 
	outfile2.write(str(num) + '\n') 
outfile2.close() 

for num in parental_fracs: 
	outfile3.write(str(num) + '\n') 
outfile3.close() 

for num in fragment_fracs: 
        outfile4.write(str(num) + '\n') 
outfile4.close() 

for num in control_fracs: 
	outfile6.write(str(num) + '\n') 
outfile6.close() 

# per-family enrichment: for each TE family, compare its per-feature percent-covered values between the two comparison pairs
# full set of families seen in any group's flanks
all_families = []
for name in group_fam_frac: 
	for feat in group_fam_frac[name]: 
		for TEid in feat: 
			if TEid not in all_families: 
				all_families.append(TEid)

comparisons = [] # each entry: [label, feature_fam_frac list for group 1, feature_fam_frac list for group 2]
comparisons.append(['parental_exons_vs_nonparental_exons', group_fam_frac['parental_exons'], group_fam_frac['nonparental_exons']])
comparisons.append(['fragments_vs_intergenic_controls', group_fam_frac['fragments'], group_fam_frac['intergenic_controls']])

fam_results = [] # each entry: [label, TEid, mean_frac_1, n_present_1, n_1, mean_frac_2, n_present_2, n_2, p_mwu, direction_mwu, p_fisher, direction_fisher]
# direction_mwu and direction_fisher are relative to group 1 (the first group named in the comparison label): enriched = more in group 1, depleted = less in group 1
mwu_pvals = []
fisher_pvals = []
for comp in comparisons: 
	label = comp[0]
	group1 = comp[1]
	group2 = comp[2]
	for TEid in all_families: 
		# per-feature percent-covered by this family for each group, with explicit 0.0 for features where the family is absent
		vec1 = []
		for feat in group1: 
			vec1.append(feat.get(TEid, 0.0))
		vec2 = []
		for feat in group2: 
			vec2.append(feat.get(TEid, 0.0))
		n1 = len(vec1)
		n2 = len(vec2)
		# mann-whitney on the percent-covered distributions (uses magnitude)
		# mannwhitneyu errors if every value is identical (eg all zeros in both groups), so guard that and call it p = 1.0
		combined = vec1 + vec2
		allsame = True
		for v in combined: 
			if v != combined[0]: 
				allsame = False
				break
		if allsame == True: 
			p_mwu = 1.0
			direction_mwu = 'equal'
		else: 
			u_mwu, p_mwu = mannwhitneyu(vec1, vec2, alternative='two-sided')
			# u_mwu is the U statistic for vec1 (group 1); U > n1*n2/2 means group 1 tends to have the larger values
			if u_mwu > (n1 * n2) / 2.0: 
				direction_mwu = 'enriched'
			elif u_mwu < (n1 * n2) / 2.0: 
				direction_mwu = 'depleted'
			else: 
				direction_mwu = 'equal'
		# fisher exact on presence/absence of this family per feature (covered > 0)
		present1 = 0
		for v in vec1: 
			if v > 0.0: 
				present1 += 1
		present2 = 0
		for v in vec2: 
			if v > 0.0: 
				present2 += 1
		oddsratio, p_fisher = fisher_exact([[present1, n1 - present1], [present2, n2 - present2]], alternative='two-sided')
		# direction from the per-feature presence rate; higher rate in group 1 = enriched in group 1
		if n1 > 0: 
			rate1 = present1 / float(n1)
		else: 
			rate1 = 0.0
		if n2 > 0: 
			rate2 = present2 / float(n2)
		else: 
			rate2 = 0.0
		if rate1 > rate2: 
			direction_fisher = 'enriched'
		elif rate1 < rate2: 
			direction_fisher = 'depleted'
		else: 
			direction_fisher = 'equal'
		fam_results.append([label, TEid, np.mean(vec1), present1, n1, np.mean(vec2), present2, n2, p_mwu, direction_mwu, p_fisher, direction_fisher])
		mwu_pvals.append(p_mwu)
		fisher_pvals.append(p_fisher)

# BH correction over all 2 * n_families tests, separately for the mwu set and the fisher set
mwu_qvals = multipletests(mwu_pvals, method='fdr_bh')[1]
fisher_qvals = multipletests(fisher_pvals, method='fdr_bh')[1]

sig_threshold = 0.05 # direction is only reported as enriched/depleted when that test's q is below this; otherwise not_significant

outfile5 = open('subtype_TE_enrichment_' + gff_type, 'w')
outfile5.write('comparison \t TE \t mean_frac_1 \t n_present_1 \t n_1 \t mean_frac_2 \t n_present_2 \t n_2 \t p_mwu \t q_mwu \t direction_mwu \t p_fisher \t q_fisher \t direction_fisher \n')
for i in range(0, len(fam_results)): 
	label = fam_results[i][0]
	TEid = fam_results[i][1]
	mean1 = fam_results[i][2]
	npres1 = fam_results[i][3]
	n1 = fam_results[i][4]
	mean2 = fam_results[i][5]
	npres2 = fam_results[i][6]
	n2 = fam_results[i][7]
	pm = fam_results[i][8]
	dm = fam_results[i][9]
	pf = fam_results[i][10]
	df = fam_results[i][11]
	qm = mwu_qvals[i]
	qf = fisher_qvals[i]
	# only report a direction if that test is significant; otherwise not_significant
	if qm < sig_threshold: 
		dm_out = dm
	else: 
		dm_out = 'not_significant'
	if qf < sig_threshold: 
		df_out = df
	else: 
		df_out = 'not_significant'
	outfile5.write(label + '\t' + TEid + '\t' + str(mean1) + '\t' + str(npres1) + '\t' + str(n1) + '\t' + str(mean2) + '\t' + str(npres2) + '\t' + str(n2) + '\t' + str(pm) + '\t' + str(qm) + '\t' + dm_out + '\t' + str(pf) + '\t' + str(qf) + '\t' + df_out + '\n')
outfile5.close()
