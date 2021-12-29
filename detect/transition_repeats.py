import numpy as np
import sys
import os

from detect.utils import genome_utils, tnf_utils, utils

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


def print_candidate_repeats(candidate_repeats, summary):
	if len(candidate_repeats) > 0:
		repeats_info = 'Misassembly detected! {} candidate repeat(s) found.\n'.format(len(candidate_repeats))
		for i,rep in enumerate(candidate_repeats):
			repeats_info += 'Candidate #{0}: ~{1}bp long, ({2}-{3}) and ({4}-{5})\n'.format(
				i+1, rep.len1, rep.start1, rep.end1, rep.start2, rep.end2)
		print(repeats_info)
		with open(summary ,'a') as fh:
			fh.write('\n\n'+repeats_info)
	else:
		print('No candidate repeats found. Stopping.')

