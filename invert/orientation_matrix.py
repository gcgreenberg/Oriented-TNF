import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

from invert.utils import genome_utils
from invert.utils import tnf_utils

def make_png(out_dir, ori_mat, chrom_len, window_len, stride):
	upper = chrom_len - (chrom_len-window_len)%stride
	extent = [0,upper,0,upper]
	extent = [(bound-window_len/2)/1e6 for bound in extent] # right edge of rightmost window
	plt.figure(figsize=(8,7))
	plt.imshow(ori_mat, origin='lower', extent=extent)
	plt.xlabel('Genome location (Mbp)', fontsize=14)
	plt.ylabel('Genome location (Mbp)',fontsize=14)
	plt.title('Orientation Decision Matrix')
	png_file = os.path.join(out_dir,'orientation_mat.png')
	plt.savefig(png_file)

def write_matrix(out_dir, ori_mat, window_locs, chrom_len):
	mat_file = os.path.join(out_dir, 'tmp', 'matrix_data.npy')
	data = np.array([ori_mat, window_locs, chrom_len], dtype=object)
	np.save(mat_file, data)

def calc_ori_mat(genome, window_len, stride):	
	tnfs,window_locs,chrom_len = tnf_utils.extract_window_tnfs(genome, window_len, stride)
	ori_mat = np.zeros((len(tnfs), len(tnfs)))
	for i,tnf1 in enumerate(tnfs):
		for j,tnf2 in enumerate(tnfs):
			ori = tnf_utils.orientation_test(tnf1,tnf2)
			ori_mat[i,j] = ori
		if i == len(tnfs)//10: print('matrix 10% complete')
		if i == len(tnfs)//2: print('matrix 50% complete')
		if i == len(tnfs)-1: print('matrix 100% complete')
	return ori_mat, window_locs, chrom_len

def make_matrix(genome_file, out_dir, chrom_id, window_len, stride, save_png, **args):
	genome = genome_utils.get_chromosome(genome_file, chrom_id=chrom_id)
	ori_mat, window_locs, chrom_len = calc_ori_mat(genome, window_len, stride)
	write_matrix(out_dir, ori_mat, window_locs, chrom_len)
	if save_png: 
		make_png(out_dir, ori_mat, chrom_len, window_len, stride)

