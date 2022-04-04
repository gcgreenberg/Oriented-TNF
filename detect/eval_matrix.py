import sys
import os
import numpy as np
from sklearn.decomposition import PCA

from detect.utils import genome_utils, tnf_utils, utils


class Range:
    def __init__(self, center, radius):
        self.center = max(center,0)
        self.lower = max(center-radius,0)
        self.upper = center+radius
    def __eq__(self, other):
        if self.center == other.center: 
            return True
        
def get_transitions(bin_pca, endpoints=False):
    trans = [i for i in range(1,len(bin_pca)) if bin_pca[i]!=bin_pca[i-1]] 
    if endpoints:
        trans.insert(0,0); trans.append(len(bin_pca))
    return trans

def binarize(ori_pca, thresh):
    bin_pca = np.sign(ori_pca - thresh)
    while len(get_transitions(bin_pca)) > 2:
        trans = get_transitions(bin_pca, endpoints=True)
        inv_lens = [trans[i]-trans[i-1] for i in range(1,len(trans))]
        blips = [i for i in range(len(inv_lens)) if inv_lens[i]==1]
        if len(blips) == 0: break
        else: bin_pca[trans[blips[0]]] *= -1
    return bin_pca

def binarize_pca(ori_pca):
    spread = max(ori_pca) - min(ori_pca)
    thresholds = np.linspace(min(ori_pca)+0.1*spread, max(ori_pca)-0.1*spread, 100)
    thresholds = thresholds[np.argsort(np.abs(thresholds))]
    for thresh in thresholds:
        bin_pca = binarize(ori_pca, thresh)
        if len(get_transitions(bin_pca)) <= 2: 
            return bin_pca
    return binarize(ori_pca, 0)

def get_trans_idx(ori_mat):
    ori_pca = PCA(n_components=1).fit_transform(ori_mat).squeeze()
    bin_pca = binarize_pca(ori_pca)
    return get_transitions(bin_pca)

def get_trans_ranges(ori_mat, chrom_len, window_len, stride):
    window_locs = utils.get_window_locs(chrom_len, window_len, stride)
    trans_idx = get_trans_idx(ori_mat)
    trans_pts = [window_locs[i]-stride//2+window_len//2 for i in trans_idx]
    trans_ranges = [Range(center=pt,radius=stride//2) for pt in trans_pts]
    return trans_ranges


def print_eval(trans_ranges, balance, summary):
    n_ranges = len(trans_ranges)
    bal_info = lambda balance: 'Genome is balanced (bal = {})'.format(balance) \
        if balance>0.8 else 'Genome is imbalanced (bal = {})'.format(balance)
    info = ''
    if n_ranges > 2:
        info += '{} transitions detected. Too many to conclude replication site locations\n'.format(n_ranges)
        if n_ranges <= 5:
            for i in range(n_ranges):
                info += 'Potential replication site #{}. Estimated location: {}\n'.format(i+1, trans_ranges[i].center)
    elif n_ranges == 2:
        info += 'Replication site #1 estimated region: {} - {}\n'.format(
            trans_ranges[0].lower, trans_ranges[0].upper)
        info += 'Replication site #2 estimated region: {} - {}\n'.format(
            trans_ranges[1].lower, trans_ranges[1].upper)
        info += 'Genome balance is {}.'.format(balance)
    elif n_ranges == 1:
        info += 'Only one transition detected. Assuming start is replication site #1\n'
        info += 'Replication site #2 estimated region: {} - {}\n'.format(
            trans_ranges[0].lower, trans_ranges[0].upper)
        info += 'Genome balance is {}'.format(balance)
    else:
        assert n_ranges == 0
        info += 'No orientation transitions detected\n'
    
    print(info)
    with open(summary, 'a') as fh:
        for i,rng in enumerate(trans_ranges): 
            fh.write('orientation transition {}: {}-{}\n'.format(i+1,rng.lower,rng.upper))

def evaluate_corrected(corr_matrix_path, stride, window_len, chrom_len, summary, **args):
    ori_mat = utils.load_file(corr_matrix_path)
    trans_ranges = get_trans_ranges(ori_mat, chrom_len, window_len, stride)
    if len(trans_ranges)==1 or len(trans_ranges)==2:
        ori, ter = utils.trans_locs(trans_ranges, chrom_len)
        balance = utils.get_balance(ori, ter, chrom_len)
    else: balance = None
    print_eval(trans_ranges, balance, summary)
            
def evaluate(matrix_path, trans_data_path, stride, window_len, chrom_len, summary, **args):
    ori_mat = utils.load_file(matrix_path)
    trans_ranges = get_trans_ranges(ori_mat, chrom_len, window_len, stride)
    if len(trans_ranges)==1 or len(trans_ranges)==2:
        ori, ter = utils.trans_locs(trans_ranges, chrom_len)
        balance = utils.get_balance(ori, ter, chrom_len)
    else: balance = None
    print_eval(trans_ranges, balance, summary)
    utils.save_file(trans_data_path, trans_ranges)
    is_imbalanced = balance<0.8 if balance is not None else False
    return is_imbalanced
