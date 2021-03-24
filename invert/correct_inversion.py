import sys
import os
import numpy as np

sys.path.append('/home/a-m/gcgreen2/code/repeat_finding/utils')
import genome_utils 
import tnf_utils

def get_best_repeat(outdir):
    repeats_file = os.path.join(outdir, 'candidate_repeats.tsv')
    try: candidate_repeats = np.loadtxt(repeats_file, dtype=int, skiprows=1)
    except: return
    num_candidates = 1 if len(np.shape(candidate_repeats))==1 else len(candidate_repeats)
    if num_candidates > 1:
        get_score = lambda rep: 10*rep[2] - abs(rep[0]-rep[4]) - abs(rep[1]-rep[5]) # 10L-dist1-dist2, dist is from expecteed to actual
        scores = [get_score(rep) for rep in candidate_repeats]
        best_rep = candidate_repeats[np.argmin(scores)]
    elif num_candidates == 1: best_rep = candidate_repeats
    return best_rep

def get_inverted_segment(chrom, repeat):
    L1 = repeat[2]; L2 = repeat[3]
    i1, i2 = np.sort([repeat[0]-1, repeat[1]-1])
    chrom_inv = chrom[:i1] + tnf_utils.rev_comp(chrom[i1:i2+L2]) + chrom[i2+L2:]
    # sanity check that start and end of repeat are the same 
    print("repeat1: {}\nrepeat2: {}".format(chrom[i1:i1+10], tnf_utils.rev_comp(chrom[i2+L2-11:i2+L2-1])))
    assert len(chrom) == len(chrom_inv)
    return chrom_inv

def get_new_chroms(file, outdir, chrom_num):
    chroms = genome_utils.get_chromosomes(file)
    repeat = get_best_repeat(outdir)
    if repeat is [] or repeat is None: return
    print('best repeat determined to be: {}'.format(repeat))
    chrom_inv = get_inverted_segment(genome_utils.get_chrom_number(chroms,chrom_num), repeat)
    new_chroms = chroms.copy()
    chrom_id = list(new_chroms)[chrom_num] 
    new_chrom_id = chrom_id + ' *** inversion-corrected ***'
    del new_chroms[chrom_id]
    new_chroms[new_chrom_id] = chrom_inv
    return new_chroms

def make_inverted_file(file, outdir, chrom_num):
    new_chroms = get_new_chroms(file, outdir, chrom_num)
    if not new_chroms:
        print('no viable repeats found. file not written')
    else:
        new_file = os.path.split(file)[1]
        prefix, suffix = os.path.splitext(new_file)
        if suffix == '.gz':
            prefix, suffix = os.path.splitext(prefix)
        if suffix == '': suffix = 'fasta'
        new_file = prefix + '_inv' + suffix
        new_file = os.path.join(outdir, new_file)
        print('making {} file with inversion-corrected chromosome at {}'.format(suffix,new_file))
        genome_utils.write_fasta(new_file, new_chroms)

def get_param(flag, default, type='int'):
    if flag in sys.argv:
        idx = sys.argv.index(flag) + 1
        param = sys.argv[idx]
    else: param = default
    if type == 'int':
        return int(param)
    elif type == 'str':
        return str(param)
    elif type == 'float':
        return float(param)

def correct    
	file = sys.argv[1]
    outdir = sys.argv[2]
    chrom_num = get_param('-chrom',0)

    make_inverted_file(file,outdir,chrom_num)






