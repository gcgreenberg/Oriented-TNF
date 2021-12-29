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

def get_trans_ranges(trans_idx, window_locs, window_len, stride, padding=0):
	trans_pts = [window_locs[i]-stride//2+window_len//2 for i in trans_idx]
	radius = stride//2 + padding
	trans_ranges = [Range(pt,radius) for pt in trans_pts]
	return trans_ranges

def project_mat(ori_mat):
	ori_pca = PCA(n_components=1).fit_transform(ori_mat).squeeze()
	ori_pca = np.sign(ori_pca)
	return ori_pca

def get_trans_idx(ori_mat):
	ori_pca = project_mat(ori_mat)
	trans_idx = [i for i in range(1,len(ori_mat)) if ori_pca[i]!=ori_pca[i-1]]
	return trans_idx

def print_transitions(trans_ranges, summary):
	n_ranges = len(trans_ranges)
	if n_ranges > 2:
		print('{} transitions detected. Too many to conclude replication site locations'.format(n_ranges))
		if n_ranges <= 5:
			for i in range(n_ranges):
				print('Potential replication site #{}. Estimated location: {}'.format(i+1, trans_ranges[i].center))
	elif n_ranges == 2:
		print('Replication site #1 estimated region: {} - {}'.format(
			trans_ranges[0].lower, trans_ranges[0].upper)) 
		print('Replication site #2 estimated region: {} - {}'.format(
			trans_ranges[1].lower, trans_ranges[1].upper))
	elif n_ranges == 1:
		print('Only one transition detected. Assuming start is replication site #1')
		print('Replication site #2 estimated region: {} - {}'.format(
			trans_ranges[0].lower, trans_ranges[0].upper))
	else:
		assert n_ranges == 0
		print('No inversions detected')
	with open(summary, 'a') as fh:
		for i,rng in enumerate(trans_ranges): 
			fh.write('inversion transition {}: {}-{}\n'.format(i+1,rng.lower,rng.upper))

def detect(matrix_data_path, trans_data_path, stride, window_len, padding, summary, **args):
	ori_mat, window_locs, seq_len = utils.load_file(matrix_data_path)
	trans_idx = get_trans_idx(ori_mat)
	trans_ranges = get_trans_ranges(trans_idx, window_locs, window_len, stride, padding)
	print_transitions(trans_ranges, summary)
	utils.save_file(trans_data_path, [trans_idx,trans_ranges])
	n_transitions = len(trans_ranges) 
	return n_transitions
