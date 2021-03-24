import os
import sys
import pandas as pd
import numpy as np
import numpy.random as rnd
import gzip
import ast
from io import StringIO
import math
from itertools import product, combinations

#UTIL_PATH = '/home/a-m/gcgreen2/code/repeat_finding/utils'
#sys.path.append(UTIL_PATH)
from invert.utils import genome_utils

BASES = {'A','T','C','G'}
BASE_PAIR = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
kmers = [''.join(kmer) for kmer in product(BASES, repeat=4)]
jmers = [''.join(kmer) for kmer in product(BASES, repeat=3)]
rev_comp = lambda seq: ''.join([BASE_PAIR[b] for b in seq[::-1]])
base_to_idx = {b:i for i,b in enumerate(BASES)}

def create_kmer_dicts():
    kmer_to_idx, idx_to_kmers = {}, {}
    cur_idx = 0
    for kmer in kmers:
        kmer_to_idx[kmer] = cur_idx
        idx_to_kmers[kmer_to_idx[kmer]] = kmer
        cur_idx += 1
    return kmer_to_idx, idx_to_kmers

def create_jmer_dicts():
    jmer_to_idx, idx_to_jmers = {}, {}
    cur_idx = 0
    for jmer in jmers:
        jmer_to_idx[jmer] = cur_idx
        idx_to_jmers[jmer_to_idx[jmer]] = jmer
        cur_idx += 1
    return jmer_to_idx, idx_to_jmers
    
kmer_to_idx, idx_to_kmers = create_kmer_dicts()
jmer_to_idx, idx_to_jmers = create_jmer_dicts()

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

def calc_rev_comps(tnfs):
    return [calc_rev_comp(tnf) for tnf in tnfs]

def get_j_prob(tnf, jmer): # j = k-1
    kmers = [jmer+base for base in BASES]
    return sum([tnf[kmer_to_idx[kmer]] for kmer in kmers])

def get_jnf(tnf): # j = k-1
    jnf = np.zeros(64)
    for jmer in jmers:
        jnf[jmer_to_idx[jmer]] = get_j_prob(tnf, jmer)
    return jnf

############## METRICS

def divergence(P, Q):
#     assert np.count_nonzero(Q) == len(Q)
    assert len(P)==len(Q)
    logadd = lambda q: 0 if q else 1e-8
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

def Euc(tnf1, tnf2):
    tot = 0
    for t1, t2 in zip(tnf1,tnf2):
        tot += (t1-t2)**2
    return np.sqrt(tot)

def orientation_test(tnf1, tnf2): # 0 if same directionality, 1 if not
    tnf_ff = (tnf1 + tnf2) / 2
    tnf_fr = (tnf1 + calc_rev_comp(tnf2)) / 2
    ent_ff = markov_entropy(tnf_ff); ent_fr = markov_entropy(tnf_fr)
    return ent_ff > ent_fr, ent_ff-ent_fr
    

########################### EXTRACTING CONTIGS

def get_contigs(genome, stride, contig_len, origin=0):
    indices = np.arange(0, len(genome)-contig_len, stride)
    contigs = []
    for idx in indices:
        contig = genome[idx: idx+contig_len]
        contigs.append(contig)
    gen_locs = np.array([int(i+origin) for i in indices])
    return contigs, gen_locs

def extract_contig_tnfs(genome, contig_len, stride, chrom_num=0):
#    contig_tnfs = []
#    chroms = genome_utils.get_chromosomes(file)
#    genome = list(chroms.values())[chrom_num]
    seq_len = len(genome)
    contigs, gen_locs = get_contigs(genome, stride, contig_len)
    contig_tnfs = [calc_tnf(contig) for contig in contigs]
    return contig_tnfs, gen_locs, seq_len
