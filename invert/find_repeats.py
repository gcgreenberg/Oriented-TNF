import numpy as np
import sys
import os
from sklearn.decomposition import PCA

sys.path.append('/home/a-m/gcgreen2/code/repeat_finding/utils')
import genome_utils
import tnf_utils

def get_dir_stem_data(dec_mat, stride, window_size=2e5):
    dec_pca = PCA(n_components=1).fit_transform(dec_mat).squeeze()
    dec_pca = np.sign(dec_pca)
    dir_pca = dec_pca.copy()
    idx_radius = int(window_size//(2*stride))
    for i in range(1,len(dir_pca)):
        cur_size = np.min([i,len(dir_pca)-i,idx_radius])
        dir_pca[i] = np.sign(np.sum(dec_pca[i-cur_size:i+cur_size])+0.1)
    return dir_pca

def get_trans_idx( dir_pca):
    N = len(dir_pca)
    trans_idx, trans_directions = [], []
    last_val = dir_pca[0]
    count = 0
    for i in range(1,N):
        count += 1
        if  last_val != dir_pca[i] and count > 3:
            trans_directions.append(dir_pca[i]-last_val)
            trans_idx.append(i)
            last_val = dir_pca[i]
            count = 0
    return trans_idx, trans_directions

def get_dec_mat_data(file, c_len, stride, chrom_num):
    contig_tnfs, gen_idx, seq_len = tnf_utils.extract_contig_tnfs(file, c_len, stride, chrom_num=chrom_num)
    dec_mat, _ = tnf_utils.calc_dir_dec_mat(contig_tnfs)
    dir_pca = get_dir_stem_data(dec_mat, stride)
    return dir_pca, seq_len

def get_repeat_regions(file, c_len, stride, chrom_num, padding=int(1e4), upper_limit=5000):
    dir_pca, seq_len = get_dec_mat_data(file, c_len, stride, chrom_num)
    N = len(dir_pca)
    if N<10: 
        print('genome too short')
        return []
    trans_idx, trans_directions = get_trans_idx(dir_pca)
    trans_pts = [int(i*stride + c_len/2 - stride/2) for i in trans_idx]
    radius = stride//2 + padding
    repeat_regions = [(pt-radius, pt+radius, pt, d) for pt,d in zip(trans_pts,trans_directions)]
    repeat_regions.insert(0,(1,1,1,0))
    repeat_regions.append((seq_len-upper_limit, seq_len, seq_len, 0))
    return repeat_regions

def parse_repeats_file(outdir, seq_id):
    repeats_file = os.path.join(outdir, 'nuc.coords')
    all_repeats = np.loadtxt(repeats_file, dtype=str, skiprows=5)
    if len(all_repeats) == 0: return []
    if len(np.shape(all_repeats))==1: all_repeats = all_repeats[np.newaxis,:]
    rc_repeats = []
    for rep in all_repeats:
        if rep[11] == seq_id and rep[12] == seq_id and rep[3] > rep[4]: # match is on same correct contig and its a rev comp match
            start1 = rep[0]; start2 = rep[4]; len1 = rep[6]; len2 = rep[7]
            rc_repeats.append([int(start1), int(start2), int(len1), int(len2)])
    return rc_repeats

def check_repeat_viable(repeat, repeat_regs):
    for i1,region1 in enumerate(repeat_regs):
        if (repeat[0] >= region1[0]) and (repeat[0] <= region1[1]):
            for i2,region2 in enumerate(repeat_regs):
                if i1 == i2: continue
               # if repeats are in regions thats correspond to opposite transitions
                elif (repeat[1] >= region2[0]) and (repeat[1] <= region2[1]) and (region1[3] != region2[3]):
                    return True, region1, region2
    return False, None, None

def get_candidate_repeats(all_repeats, repeat_regs):
    candidate_repeats = []
    if len(repeat_regs) > 6 or len(repeat_regs) <= 4 : return candidate_repeats # too many regions indicates a messy genome, too little means a normal looking one
    for repeat in all_repeats:
        repeat_viable, region1, region2 = check_repeat_viable(repeat, repeat_regs)
        if repeat_viable:    
         # start1      start2     len1       len2       expected1    expected2
            full_repeat = (repeat[0], repeat[1], repeat[2], repeat[3], region1[2], region2[2])
            candidate_repeats.append(full_repeat)
    return candidate_repeats

def make_repeat_regs_file(repeat_regs, outdir):
    regs_file = os.path.join(outdir, 'repeat_regions.tsv')
    with open(regs_file, 'w') as fh:
        fh.write("lower\tupper\texpected\ttransition_type\n")
        for reg in repeat_regs:
            fh.write("{}\t{}\t{}\t{}\n".format(*reg))

def make_candidate_repeats_file(file, outdir, c_len, stride, chrom_num, padding, seq_id):
    repeat_regs = get_repeat_regions(file, c_len, stride, chrom_num, padding=padding)
    print('potential regions for repeats: {}'.format(repeat_regs))
    make_repeat_regs_file(repeat_regs, outdir)
    all_repeats = parse_repeats_file(outdir, seq_id)    
    candidate_repeats = get_candidate_repeats(all_repeats, repeat_regs)
    if len(candidate_repeats) > 0:
        repeats_file = os.path.join(outdir, 'candidate_repeats.tsv')
        with open(repeats_file, 'w') as fh:
            fh.write("seq_f\tseq_r\tlength1\tlength2\texpected_f\texpected_r\n")
            for rep in candidate_repeats:
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(*rep))
        with open(repeats_file, 'r') as fh: 
            print(fh.read()+"\n")
            print('candidate repeats written to {}'.format(repeats_file))
    else:
        print('no viable repeats found. stopping.')
 
def print_info(file, c_len, stride, outdir, chrom_num):
    chroms = genome_utils.get_chromosomes(file)
    print('file: {}'.format(file))
    print('gen len: {}, num chrom: {}'.format(len(genome_utils.get_chrom_number(chroms,0)), len(chroms)))
    print('c_len: {},  stride: {}'.format(c_len, stride))
    print('out directory: {}'.format(outdir))
    seq_id = list(chroms.keys())[chrom_num]
    seq_id = seq_id.split(' ')[0]
    print('seq_id: {}'.format(seq_id))
    return seq_id

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

if __name__ == '__main__':
    file = sys.argv[1]
    outdir = sys.argv[2]
    c_len = get_param('-len', 5e4)
    stride = get_param('-stride', 4e4)
    chrom_num = get_param('-chrom',0)
    padding = get_param('-pad',1e4)

    seq_id = print_info(file, c_len, stride, outdir, chrom_num)
    make_candidate_repeats_file(file,outdir,c_len,stride,chrom_num,padding,seq_id) 
