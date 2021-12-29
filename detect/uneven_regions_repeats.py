import numpy as np
import sys
import os

from detect.utils import genome_utils, tnf_utils, utils
 
    
def is_candidate(rep, lower, upper, chrom_len):
    rep_in_ranges = (repeat.start1 > lower and repeat.end1 < upper and (repeat.end2 > upper or repeat.start2 < lower)) or \
        (repeat.start2 > lower and repeat.end2 < upper and (repeat.end1 > upper or repeat.start1 < lower))
    if rep_in_ranges:
        old_reg_diff = abs(chrom_len - 2*abs(upper - lower))
        new_reg_diff = abs(chrom_len - 2*abs(upper + lower - rep.start1 - rep.start2))
        return new_reg_diff < old_reg_diff
    return False

def print_candidate_repeats(candidate_repeats, summary):
	if len(candidate_repeats) > 0:
		repeats_info = 'Misassembly detected! {} candidate repeat found.\n'.format(len(candidate_repeats))
		for i,rep in enumerate(candidate_repeats):
			repeats_info += 'Candidate #{0}: ~{1}bp long, ({2}-{3}) and ({4}-{5})\n'.format(
				i+1, rep.len1, rep.start1, rep.end1, rep.start2, rep.end2)
	else:
		repeats_info = 'No candidate repeats found. Stopping.'
	print(repeats_info)
	with open(summary ,'a') as fh:
		fh.write('\n\n'+repeats_info)
        
def get_inversion_locs(trans_ranges, chrom_len):
    lower = trans_ranges[0].center
    upper = trans_ranges[1].center if len(trans_ranges)==2 else chrom_len
    return lower, upper
    

def find_candidate_repeats(all_repeats, trans_data_path, matrix_data_path, candidates_path, **args):
    _,trans_ranges = utils.load_file(trans_data_path)
    _,_,chrom_len = utils.load_file(matrix_data_path)
    lower, upper, reg_diff = get_inversion_info(trans_ranges, chrom_len)
    candidate_reps = [rep for rep in all_repeats if is_candidate(rep, lower, upper, chrom_len)]
    if len(candidate_reps) > 0: utils.save_file(candidates_path, candidate_reps)
    