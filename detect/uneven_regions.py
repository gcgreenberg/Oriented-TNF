import numpy as np
import sys
import os

from detect.utils import genome_utils, tnf_utils, utils

def is_candidate(rep, trans_ranges, chrom_len, len_thresh): # ori < ter
    ori, ter = utils.trans_locs(trans_ranges, chrom_len)
    new_balance = utils.new_balance(rep, ori, ter, chrom_len) 
    sufficient_len = rep.len1 > len_thresh
    return new_balance>0.8 and sufficient_len

def print_candidate_repeats(candidate_repeats, summary):
    if len(candidate_repeats) > 0:
        repeats_info = 'Misassembly detected! {} candidate repeat found.\n'.format(len(candidate_repeats))
        for i,rep in enumerate(candidate_repeats):
            repeats_info += 'Candidate #{0}: ~{1}bp long, locations ({2}-{3}) and ({4}-{5})\n'.format(
                i+1, rep.len1, rep.start1, rep.end1, rep.start2, rep.end2)
    else:
        repeats_info = 'No candidate repeats found. Stopping.'
    print(repeats_info)
    with open(summary ,'a') as fh:
        fh.write('\n\n'+repeats_info)
    
def find_candidate_repeats(all_repeats, trans_data_path, candidates_path, chrom_len, len_thresh, summary, **args):
    trans_ranges = utils.load_file(trans_data_path)
    candidate_reps = [rep for rep in all_repeats if is_candidate(rep, trans_ranges, chrom_len, len_thresh)]
    print_candidate_repeats(candidate_reps, summary)
    if len(candidate_reps) > 0: utils.save_file(candidates_path, candidate_reps)
    return len(candidate_reps) > 0
    