import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from detect.utils import genome_utils, tnf_utils, utils

def make_png(png_path, ori_mat, chrom_len):
	extent = [0,chrom_len,0,chrom_len]
	plt.figure(figsize=(8,7))
	plt.imshow(ori_mat, origin='lower', extent=extent)
	plt.xlabel('Genome location (bp)', fontsize=14)
	plt.ylabel('Genome location (bp)',fontsize=14)
	plt.title('Orientation Decision Heatmap')
	plt.savefig(png_path)

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

def make_matrix(genome_path, chrom_id, window_len, stride, png_path, matrix_data_path, **args):
	genome = genome_utils.get_chromosome(genome_path, chrom_id=chrom_id)
	ori_mat, window_locs, chrom_len = calc_ori_mat(genome, window_len, stride)
	utils.save_file(matrix_data_path, [ori_mat, window_locs, chrom_len])
	make_png(png_path, ori_mat, chrom_len)