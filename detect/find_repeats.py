import numpy as np
import sys
import os

from detect.utils import genome_utils, tnf_utils, utils

class Repeat:
	def __init__(self, rep):
		self.start1 = int(rep[0])
		self.start2 = int(rep[3])
		self.len1 = int(rep[6])
		self.len2 = int(rep[7])
		self.end1 = int(rep[1])
		self.end2 = int(rep[4])
		self.idy = float(rep[9])
		self.seq_id1 = rep[17]
		self.seq_id2 = rep[18]

def parse_repeats_file(repeats_path):
	all_repeats = np.loadtxt(repeats_path, dtype=str, skiprows=5)
	if len(all_repeats) == 0: return
	if len(np.shape(all_repeats))==1: all_repeats = all_repeats[np.newaxis,:]
	all_repeats = [Repeat(rep) for rep in all_repeats]
	return all_repeats

def is_repeat_viable(repeat, chrom_id):
	return repeat.seq_id1 == chrom_id and \
		repeat.seq_id2 == chrom_id and \
		repeat.end2 < repeat.start2 and \
		repeat.idy >= 95

def rep_in_ranges(repeat, trans_ranges):
	for i1,range1 in enumerate(trans_ranges):
		if (repeat.start1 >= range1.lower) and (repeat.start1 <= range1.upper):
			for i2,range2 in enumerate(trans_ranges):
				if i1 == i2: continue
			   # if repeats are in regions thats correspond to opposite transitions
				if (repeat.start2 >= range2.lower) and \
						(repeat.start2 <= range2.upper) and \
						(range1 != range2):
					return True
	return False

def remove_duplicates(repeats):
	thresh = 1000
	is_duplicate = lambda rep1, rep2: (abs(rep1.start1-rep2.start1) < thresh and \
		abs(rep1.start2-rep2.start2) < thresh) or \
		(abs(rep1.start1-rep2.end2) < thresh and \
		abs(rep1.start2-rep2.end1) < thresh)
	candidate_repeats = []
	for i,rep1 in enumerate(repeats):
		if np.all([not is_duplicate(rep1, rep2) for rep2 in candidate_repeats[:i]]):
			candidate_repeats.append(rep1)
	return candidate_repeats

def get_candidate_repeats(repeats_path, trans_data_path, chrom_id):
	trans_ranges = utils.load_file(trans_data_path)
	all_repeats = parse_repeats_file(repeats_path)
	candidate_repeats = [repeat for repeat in all_repeats if \
			is_repeat_viable(repeat, chrom_id) and rep_in_ranges(repeat, trans_ranges)]
	candidate_repeats = remove_duplicates(candidate_repeats)
	return candidate_repeats

def print_candidate_repeats(candidate_repeats):
	if len(candidate_repeats) > 0:
		print('Misassembly detected! {} candidate repeat(s) found.'.format(
			len(candidate_repeats)))
		for i,rep in enumerate(candidate_repeats):
			print('Candidate #{0}: ~{1}bp long, ({2}-{3}) and ({4}-{5})'.format(
				i, rep.len1, rep.start1, rep.end1, rep.start2, rep.end2))
	else:
		print('No candidate repeats found. Stopping.')

def find_candidate_repeats(repeats_path, chrom_id, trans_data_path, **args):
	candidate_repeats = get_candidate_repeats(repeats_path, trans_data_path, chrom_id)
	print_candidate_repeats(candidate_repeats)
	utils.save_file(args['candidates_path'], candidate_repeats)
	return len(candidate_repeats) > 0

def run_nucmer(mummer_path, delta_path, repeats_path, tmp_genome_path, **args):
	nuc_path = os.path.join(mummer_path, 'nucmer')
	show_coords_path = os.path.join(mummer_path, 'show-coords')
	os.system('{0} -r -l 15 -c 50 -g 20 -b 100 --delta {1} {2} {2}'.format(nuc_path, delta_path, tmp_genome_path))
	os.system('{0} -r -l -c {1} > {2}'.format(show_coords_path, delta_path, repeats_path))
