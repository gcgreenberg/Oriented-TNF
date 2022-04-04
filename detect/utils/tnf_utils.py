import os
import sys
import numpy as np
import numpy.random as rnd
from itertools import product

from detect.utils import genome_utils

BASES = {'A','T','C','G'}
BASE_PAIR = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
kmers = [''.join(kmer) for kmer in product(BASES, repeat=4)]
jmers = [''.join(kmer) for kmer in product(BASES, repeat=3)]
rev_comp = lambda seq: ''.join([BASE_PAIR[b] for b in seq[::-1]])
base_to_idx = {b:i for i,b in enumerate(BASES)}

def create_dicts(kmers):
    kmer_to_idx, idx_to_kmers = {}, {}
    cur_idx = 0
    for kmer in kmers:
        kmer_to_idx[kmer] = cur_idx
        idx_to_kmers[kmer_to_idx[kmer]] = kmer
        cur_idx += 1
    return kmer_to_idx, idx_to_kmers
 
kmer_to_idx, idx_to_kmers = create_dicts(kmers)
jmer_to_idx, idx_to_jmers = create_dicts(jmers)

def calc_tnf(seq):
    tnf = np.zeros(256)
    for i in range(len(seq)-4+1):
        kmer = seq[i:i+4]
        tnf[kmer_to_idx[kmer]] += 1
    return tnf / np.sum(tnf)
    
def calc_rev_comp(tnf):
    rc_tnf = np.zeros_like(tnf)
    for p,kmer in zip(tnf,kmers):
        idx = kmers.index(rev_comp(kmer))
        rc_tnf[idx] = p
    return rc_tnf

def get_jnf(tnf): # j = k-1
    jnf = np.zeros(64)
    for jmer in jmers:
        kmers = [jmer+base for base in BASES]
        jmer_prob = sum([tnf[kmer_to_idx[kmer]] for kmer in kmers])
        jnf[jmer_to_idx[jmer]] = jmer_prob
    return jnf

############## METRICS

def divergence(P, Q):
    assert len(P)==len(Q)
    logadd = lambda q: 0 if q else 1e-12
    return sum([p * np.log2((p+1e-12)/(q+logadd(q))) for p,q in zip(P,Q)]) 

def entropy(P):
    assert np.isclose(sum(P), 1)
    logadd = lambda p: 0 if p else min(P[P!=0])/2
    return -sum([p * np.log2((p+1e-12)) for p in P]) 
    
def markov_divergence(Pk,Qk):
    Pj, Qj = get_jnf(Pk), get_jnf(Qk)
    return divergence(Pk,Qk) - divergence(Pj,Qj)

def markov_entropy(Pk):
    Pj = get_jnf(Pk)
    return entropy(Pk) - entropy(Pj)

def orientation_test(tnf1, tnf2): # 0 if same orientation, 1 if not
    tnf_same = (tnf1 + tnf2) / 2
    tnf_opposite = (tnf1 + calc_rev_comp(tnf2)) / 2
    ent_same = markov_entropy(tnf_same)
    ent_opposite = markov_entropy(tnf_opposite)
    return ent_same > ent_opposite
    
########################### EXTRACTING CONTIGS

def extract_window_tnfs(genome, window_len, stride, chrom_num=0):
    gen_locs = np.arange(0, len(genome)-window_len, stride)
    window_tnfs = []
    for idx in gen_locs:
        window = genome[idx: idx+window_len]
        window_tnfs.append(calc_tnf(window))
    return window_tnfs, gen_locs
