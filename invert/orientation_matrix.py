import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from sklearn.decomposition import PCA

from invert.utils import genome_utils
from invert.utils import tnf_utils


def get_tick_len(seq_len):
    if seq_len > 3e6: tick_len=5e5
    elif seq_len > 1e6: tick_len=1e5
    elif seq_len > 3e5: tick_len = 5e4
    elif seq_len > 3e4: tick_len = 5e3
    else: tick_len = 5e2
    return int(tick_len),int(np.floor(np.log10(tick_len)))

def get_heatmap_axis_labels(gen_idxs, seq_lens):
    plot_gen_idxs, tnf_idxs = [],[]
    for gen_idx,seq_len in zip(gen_idxs, seq_lens):
        tick_len,tick_exp = get_tick_len(seq_len)
        plot_gen_idx = np.arange(round(gen_idx[0],-tick_exp),gen_idx[-1],tick_len)
        tnf_idx = [find_interpolated_idx(gen_idx, i) for i in plot_gen_idx]
        plot_gen_idxs.append(plot_gen_idx), tnf_idxs.append(tnf_idx)
    return plot_gen_idxs, tnf_idxs




def plot_comp_dec_mat(D, chrom_len, window_len=10000, stride=5000):
#    plot_gen_idxs, tnf_idxs = get_heatmap_axis_labels(gen_idxs,seq_lens)
	plt.figure(figsize=(8,7))
	extent = (np.array([0,chrom_len,0,chrom_len])-window_len//2) / 1e6
	plt.imshow(D, origin='lower',extent=extent)
	#plt.xticks(tnf_idxs[1], plot_gen_idxs[1],rotation=90); plt.yticks(tnf_idxs[0], plot_gen_idxs[0])
	plt.xlabel('Genome location (Mbp)',fontsize=14); plt.ylabel('Genome location (Mbp)',fontsize=14)
	plt.title('Orientation Decision Matrix')


def make_png(D,out_dir,chrom_len,window_len,stride):
	png_file = os.path.join(out_dir,'orientation_mat.png')
	plot_comp_dec_mat(D,chrom_len,window_len,stride)
	plt.savefig(png_file)


def write_matrix(orientation_mat, out_dir):
	mat_file = os.path.join(out_dir, 'orientation_mat.npy')
	np.save(mat_file, orientation_mat)


def calc_ori_mat(genome, window_len, stride):	
	tnfs,_,chrom_len = tnf_utils.extract_contig_tnfs(genome, window_len, stride)
	D = np.zeros((len(tnfs), len(tnfs)))
	score_mat = np.zeros_like(D)
	for i,tnf1 in enumerate(tnfs):
		if i == len(tnfs)//10: print('matrix 10% complete')
		if i == len(tnfs)//2: print('matrix 50% complete')
		for j,tnf2 in enumerate(tnfs):
			dir, score = tnf_utils.orientation_test(tnf1,tnf2)
			D[i,j] = dir
			score_mat[i,j] = score
	return D, chrom_len


def make_matrix(genome_file, out_dir, **args):
	genome = genome_utils.get_chromosome(genome_file, chrom_id=args['chrom_id'])
	orientation_mat,chrom_len = calc_ori_mat(genome, window_len=args['window_len'], stride=args['stride'])
	write_matrix(orientation_mat, out_dir)
	if args['save_png']: make_png(orientation_mat,out_dir,chrom_len,window_len=args['window_len'],stride=args['stride'])

