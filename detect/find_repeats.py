import numpy as np
import sys
import os

from detect.utils import genome_utils, tnf_utils, utils

class Repeat:
    def __init__(self, rep):
        self.start1 = int(rep[0])
        self.start2 = int(rep[2])
        self.len1 = int(rep[4])
        self.len2 = int(rep[5])
        self.end1 = int(rep[1])
        self.end2 = int(rep[3])
        self.idy = float(rep[6])
        self.seq_id1 = rep[7]
        self.seq_id2 = rep[8]

def parse_repeats_file(repeats_path):
    n_lines = sum(1 for line in open(repeats_path))
    if n_lines == 4: return [] # header is 4 lines
    all_repeats = np.loadtxt(repeats_path, delimiter='\t', dtype=str, skiprows=4)
    if len(np.shape(all_repeats))==1: all_repeats = all_repeats[np.newaxis,:]
    all_repeats = [Repeat(rep) for rep in all_repeats]
    return all_repeats

def is_repeat_viable(repeat, chrom_id):
    return repeat.seq_id1 == chrom_id and \
        repeat.seq_id2 == chrom_id and \
        repeat.end2 < repeat.start2

def remove_duplicates(repeats):
    thresh = 100
    is_duplicate = lambda rep1, rep2: \
        (abs(rep1.start1-rep2.start1)<thresh and abs(rep1.start2-rep2.start2)<thresh and abs(rep1.len1-rep2.len1)<thresh) or \
        (abs(rep1.start1-rep2.end2)<thresh and abs(rep1.start2-rep2.end1)<thresh and abs(rep1.len1-rep2.len1)<thresh)
    all_repeats = []
    for i,rep1 in enumerate(repeats):
        if np.all([not is_duplicate(rep1, rep2) for rep2 in all_repeats[:i]]):
            all_repeats.append(rep1)
    return all_repeats


def find_all_repeats(repeats_path, chrom_id, summary, **args):
    all_repeats = parse_repeats_file(repeats_path)
    all_repeats = [repeat for repeat in all_repeats if is_repeat_viable(repeat, chrom_id)]
    all_repeats = remove_duplicates(all_repeats)
    return all_repeats

def run_nucmer(mummer_path, delta_path, repeats_path, tmp_genome_path, **args):
    nuc_path = os.path.join(mummer_path, 'nucmer')
    show_coords_path = os.path.join(mummer_path, 'show-coords')
    os.system('{0} -r --delta {1} {2} {2}'.format(nuc_path, delta_path, tmp_genome_path))
    os.system('{0} -r -T -H -I 95 {1} > {2}'.format(show_coords_path, delta_path, repeats_path))
