# usage: python reorient_reads.py (infasta) (outfastaname)
import sys 

outfasta_name = sys.argv[2]

def read_fasta(path):
	"""Yield (header, sequence) tuples from a fasta file."""
	header = None
	seq_parts = []
	with open(path) as f:
		for line in f:
			line = line.rstrip('\n')
			if line.startswith('>'):
				if header is not None:
					yield header, ''.join(seq_parts)
				header = line[1:]
				seq_parts = []
			else:
				seq_parts.append(line)
		if header is not None:
			yield header, ''.join(seq_parts)

def longest_homopolymer_run(seq, base):
	"""Return the length of the longest consecutive run of `base` in seq."""
	seq = seq.upper()
	base = base.upper()
	longest = 0
	current = 0
	for c in seq:
		if c == base:
			current += 1
			if current > longest:
				longest = current
		else:
			current = 0
	return longest

def revcomp(seq):
	"""Return the reverse complement of a DNA sequence."""
	complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
	return seq.translate(complement)[::-1]

def write_fasta_record(out, header, seq, width=60):
	"""Write one fasta record with sequence wrapped to `width` chars per line."""
	out.write('>' + header + '\n')
	for i in range(0, len(seq), width):
		out.write(seq[i:i+width] + '\n')

outfile2 = open(sys.argv[1] + '_TS_STRANDS', 'w')

with open(outfasta_name, 'w') as out:
	for header, seq in read_fasta(sys.argv[1]):
		head = seq[:200]
		tail = seq[-200:]
		t_run_5 = longest_homopolymer_run(head, 'T')
		a_run_3 = longest_homopolymer_run(tail, 'A')
		if t_run_5 > a_run_3: 
			newseq = revcomp(seq) 
			outfile2.write(header + '\t' + '-' + '\n')
		else: 
			newseq = seq
			outfile2.write(header + '\t' + '+' + '\n')
		write_fasta_record(out, header, newseq)
outfile2.close()
