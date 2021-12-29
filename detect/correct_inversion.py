import sys
import os
import numpy as np

from detect import orientation_matrix
from detect.utils import genome_utils, tnf_utils, utils

def get_inverted_chrom(chrom, repeat):
	_,i1,i2,_ = np.sort([repeat.start1, repeat.start2, repeat.end1, repeat.end2])
	i2 -= 1
	chrom_inv = chrom[:i1] + tnf_utils.rev_comp(chrom[i1:i2]) + chrom[i2:]
	assert len(chrom) == len(chrom_inv)
	return chrom_inv, i1, i2

def correct_misassembly(old_chroms, repeat, chrom_id):
	chrom_inv, i1, i2 = get_inverted_chrom(old_chroms[chrom_id], repeat)
	corrected_chroms = old_chroms.copy()
	new_chrom_id = chrom_id + ' *** inversion-corrected ({}-{}) ***'.format(i1,i2)
	del corrected_chroms[chrom_id]
	corrected_chroms[new_chrom_id] = chrom_inv
	return corrected_chroms, chrom_inv

def get_longest_repeat(candidate_repeats):
	lens = [r.len1 for r in candidate_repeats]
	return candidate_repeats[np.argmax(lens)]

def save_fasta(chroms, **args):
	corrected_genome_path = utils.get_corrected_genome_path(**args)
	print('making inversion-corrected chromosome, filepath: {}'.format(corrected_genome_path))
	genome_utils.write_fasta(corrected_genome_path, chroms)

def orientation_mat_png(chrom_inv, window_len, stride, **args):
	corrected_png_path = utils.get_corrected_png_path(**args)
	ori_mat, _, chrom_len = orientation_matrix.calc_ori_mat(chrom_inv, window_len, stride)
	orientation_matrix.make_png(corrected_png_path, ori_mat, chrom_len)

def correct(tmp_genome_path, candidates_path, chrom_id, window_len, stride, **args):
	old_chroms = genome_utils.get_chromosomes(tmp_genome_path)
	candidate_repeats = utils.load_file(candidates_path)
	best_repeat = get_longest_repeat(candidate_repeats)
	utils.print_banner('CORRECTING MISASSEMBLY')
	corrected_chroms, chrom_inv = correct_misassembly(old_chroms, best_repeat, chrom_id)
	save_fasta(chroms=corrected_chroms, **args)
	utils.print_banner('MAKING CORRECTED MATRIX')
	orientation_mat_png(chrom_inv, window_len, stride, **args)
