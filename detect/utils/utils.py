import sys
import os
import numpy as np

from detect.utils import genome_utils

def print_banner(x):
    print('\n=================== {} ==================='.format(str(x)))

def print_header():
    print_banner('ORIENTED - TNF')
    for _ in range(3): print('=======================================')
        
def save_file(filepath, x):
    x = np.array(x, dtype=object)
    np.save(filepath, x, allow_pickle=True)
    
def load_file(filepath):
    return np.load(filepath, allow_pickle=True)

def print_matrix_info(genome_path, window_len, stride, out_dir, chrom_id, summary, **args):
    info = 'genome file: {}\n'.format(genome_path) + \
            'sequence id chosen: {}\n'.format(chrom_id) + \
            'window length: {},  stride: {}\n'.format(window_len, stride) + \
            'output directory: {}\n'.format(out_dir)
    print(info)
    with open(summary, 'a') as fh:
        fh.write(info)
    
def setup_out_dir(out_dir, genome_path, tmp_genome_path, **args):
    tmp_dir = os.path.join(out_dir, 'tmp')
    os.makedirs(tmp_dir, exist_ok=True)
    chroms = genome_utils.get_chromosomes(genome_path)
    genome_utils.write_fasta(tmp_genome_path, chroms)
    
def get_corrected_genome_path(genome_path, out_dir, **args):
#     index = '' if index==0 else str(index+1)
    corrected_path = os.path.split(genome_path)[1]
    prefix, suffix = os.path.splitext(corrected_path)
    if suffix == '.gz':
        prefix, suffix = os.path.splitext(prefix)
    if suffix == '': suffix = 'fasta'
    corrected_path = prefix + '_corrected' + suffix
    return os.path.join(out_dir, corrected_path)

def trans_locs(trans_ranges, chrom_len):
    t1 = trans_ranges[0].center
    t2 = trans_ranges[1].center if len(trans_ranges)==2 else chrom_len
    ori, ter = np.sort([t1,t2])
    return ori, ter

def get_balance(ori, ter, chrom_len):
    balance = 1-abs(chrom_len-2*(ter-ori))/chrom_len
    return round(balance, 3)

def new_balance(rep, ori, ter, chrom_len):
    start1, end1, start2, end2 = np.sort([rep.start1, rep.end1, rep.start2, rep.end2]) 
    balance = -1
    if end1 < ori and start2 > ori and end2 < ter:
        new_ori = end1+(start2-ori); new_ter = ter
        balance = get_balance(new_ori, new_ter, chrom_len)
    elif start1 > ori and end1 < ter and start2 > ter:
        new_ori = ori; new_ter = end1+(start2-ter)
        balance = get_balance(new_ori, new_ter, chrom_len)
    return balance

def get_window_locs(chrom_len, window_len, stride):
    return np.arange(0, chrom_len-window_len, stride)
