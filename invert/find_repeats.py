import numpy as np
import sys
import os
from sklearn.decomposition import PCA

from invert.utils import genome_utils
from invert.utils import tnf_utils

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

def parse_repeats_file(out_dir, chrom_id):
	repeats_file = os.path.join(out_dir, 'tmp', 'nuc.coords')
	all_repeats = np.loadtxt(repeats_file, dtype=str, skiprows=5)
	if len(all_repeats) == 0: return
	if len(np.shape(all_repeats))==1: all_repeats = all_repeats[np.newaxis,:]
	all_repeats = [Repeat(rep) for rep in all_repeats]
	return all_repeats

def is_repeat_viable(repeat, chrom_id):
	if repeat.seq_id1 == chrom_id and \
			repeat.seq_id2 == chrom_id and \
			repeat.end2 < repeat.start2 and \
			repeat.idy >= 95:
			return True
	else: return False

def rep_in_ranges(repeat, trans_ranges):
	for i1,range1 in enumerate(trans_ranges):
		if (repeat.start1 >= range1.lower) and (repeat.start1 <= range1.upper):
			for i2,range2 in enumerate(trans_ranges):
				if i1 == i2: continue
				if repeat.start1==2630454: print
			   # if repeats are in regions thats correspond to opposite transitions
				if (repeat.start2 >= range2.lower) and \
						(repeat.start2 <= range2.upper) and \
						(range1 != range2):
					return True
	return False

def get_best_candidate(candidate_repeats):
	lens = [rep.len1 for rep in candidate_repeats]
	return candidate_repeats[np.argmax(lens)]

def get_candidate_repeats(all_repeats, trans_ranges, chrom_id):
	candidate_repeats = [repeat for repeat in all_repeats if \
			is_repeat_viable(repeat, chrom_id) and rep_in_ranges(repeat, trans_ranges)]
#	print(chrom_id)
#	for rep in all_repeats: print('{}, {}, {}, {}'.format(rep.seq_id1,rep.seq_id2,rep.end2 < rep.start2,rep.idy))
	#print([rep for rep in all_repeats if is_repeat_viable(rep,chrom_id)])
	if len(candidate_repeats) > 0: return candidate_repeats

def print_candidate_repeats(candidate_repeats, best_repeat):
	print('Misassembly detected! {} candidate repeat(s) found.'.format(len(candidate_repeats)))
	print('Best repeat is {}bp long, and starts at locations {} and {}'.format(
		best_repeat.len1, best_repeat.start1, best_repeat.start2))

def find_best(trans_ranges, out_dir, chrom_id, **args):
	all_repeats = parse_repeats_file(out_dir, chrom_id)
	candidate_repeats = get_candidate_repeats(all_repeats, trans_ranges, chrom_id)
	if candidate_repeats is not None:
		best_repeat = get_best_candidate(candidate_repeats)
		print_candidate_repeats(candidate_repeats, best_repeat)
		return best_repeat
	else:
		print('no candidate repeats found. stopping.')
