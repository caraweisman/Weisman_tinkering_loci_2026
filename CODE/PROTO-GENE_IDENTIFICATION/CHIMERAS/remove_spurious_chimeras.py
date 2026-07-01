import sys
import numpy as np

# usage: remove_spurious_chimeras.py (blastout info from orfs)

blastout = np.genfromtxt(sys.argv[1], dtype=str, delimiter='\t')

orf_blastout_dict = {}
for row in blastout:
	orf_blastout_dict.setdefault(row[0], []).append(row)


def overlap(iv1, iv2):
	a, b = iv1
	c, d = iv2
	return max(0, min(b, d) - max(a, c))


def filter_contained_indices(intervals, max_aa=50, frac=0.8):
	"""Drop interval i if some bigger interval overlaps it by > max_aa OR > frac of its length."""
	def is_swallowed(i):
		a, b = intervals[i]
		length = b - a
		if length == 0:
			return False
		for j, (c, d) in enumerate(intervals):
			if i == j:
				continue
			if (d - c) <= length:
				continue
			ov = overlap(intervals[i], intervals[j])
			if ov > max_aa or ov / length > frac:
				return True
		return False
	return [i for i in range(len(intervals)) if not is_swallowed(i)]


def dedupe_by_interval(rows):
	"""Keep the row with the best (lowest) e-value per (hitstart, hitend)."""
	best = {}
	for r in rows:
		key = (int(r[2]), int(r[3]))
		evalue = float(r[8])
		if key not in best or evalue < float(best[key][8]):
			best[key] = r
	return list(best.values())


with open(sys.argv[1] + '_filtered', 'w') as outfile:
	for rows in orf_blastout_dict.values():
		rows = dedupe_by_interval(rows)
		hit_intervals = [(int(r[2]), int(r[3])) for r in rows]
		keep = filter_contained_indices(hit_intervals, max_aa=50, frac=0.8)
		if len(keep) <= 1:
			continue
		for i in keep:
			outfile.write('\t'.join(rows[i]) + '\n')
