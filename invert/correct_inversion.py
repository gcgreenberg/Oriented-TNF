import sys
import os
import numpy as np

from invert.utils import genome_utils 
from invert.utils import tnf_utils

def get_inverted_chrom(chrom, repeat):
	_,i1,i2,_ = np.sort([repeat.start1, repeat.start2, repeat.end1, repeat.end2])
	i2 -= 1
	chrom_inv = chrom[:i1] + tnf_utils.rev_comp(chrom[i1:i2]) + chrom[i2:]
	# sanity check that start and end of repeat are the same 
	assert len(chrom) == len(chrom_inv)
	return chrom_inv

def reinvert_misassembly(old_chroms, best_repeat, chrom_id):
	chrom_inv = get_inverted_chrom(old_chroms[chrom_id], best_repeat)
	corrected_chroms = old_chroms.copy()
	new_chrom_id = chrom_id + ' *** inversion-corrected ***'
	del corrected_chroms[chrom_id]
	corrected_chroms[new_chrom_id] = chrom_inv
	return corrected_chroms

def write_corrected_fasta(genome_file, outdir, corrected_chroms):
	new_file = os.path.split(genome_file)[1]
	prefix, suffix = os.path.splitext(new_file)
	if suffix == '.gz':
		prefix, suffix = os.path.splitext(prefix)
	if suffix == '': suffix = 'fasta'
	new_file = prefix + '_inv' + suffix
	new_file = os.path.join(outdir, new_file)
	print('making inversion-corrected chromosome, filename {}'.format(new_file))
	genome_utils.write_fasta(new_file, corrected_chroms)
	return new_file

def correct(genome_file, out_dir, best_repeat, chrom_id, **args):
	old_chroms = genome_utils.get_chromosomes(genome_file)
	corrected_chroms = reinvert_misassembly(old_chroms, best_repeat, chrom_id)
	new_file = write_corrected_fasta(genome_file, out_dir, corrected_chroms)
	return new_file

