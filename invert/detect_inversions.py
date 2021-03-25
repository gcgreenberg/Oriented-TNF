import sys
import os
import numpy as np
from sklearn.decomposition import PCA

sys.path.append('/home/a-m/gcgreen2/code/repeat_finding/utils')
from invert.utils import genome_utils
from invert.utils import tnf_utils


class Range:
	def __init__(self, center, radius):
		self.center = max(center,0)
		self.lower = max(center-radius,0)
		self.upper = center+radius

def get_trans_ranges(trans_idx, window_locs, window_len, stride, padding=0):
	trans_pts = [window_locs[i]-stride//2+window_len//2 for i in trans_idx]
	radius = stride//2 + padding
	trans_ranges = [Range(pt,radius) for pt in trans_pts]
	print(trans_idx,trans_pts,radius,trans_ranges)
	return trans_ranges

def project_mat(ori_mat):
	ori_pca = PCA(n_components=1).fit_transform(ori_mat).squeeze()
	ori_pca = np.sign(ori_pca)
	return ori_pca

def get_trans_idx(ori_mat):
	ori_pca = project_mat(ori_mat)
	trans_idx = [i for i in range(1,len(ori_mat)) if ori_pca[i]!=ori_pca[i-1]]
	return trans_idx

def print_transitions(trans_ranges, seq_len):
	n_ranges = len(trans_ranges)
	if n_ranges > 5:
		print('Too many inversion transitions detected')
	elif n_ranges > 2:
		print('{} transitions detected. Too many to find replication sites'.format(n_ranges))
		for i in range(n_ranges//2):
			print('Inversion {} estimated location: {} - {}'.format(
				i+1, trans_ranges[2*i].center, trans_ranges[2*i+1].center))
		if n_ranges%2 == 1:
			print('Inversion {} estimated location: {} - end ({})'.format(
				(n_ranges+1)//2, trans_ranges[-1].center, seq_len))	
	elif n_ranges == 2:
		print('Replication site 1 estimated region: {} - {}'.format(
			trans_ranges[0].lower, trans_ranges[0].upper)) 
		print('Replication site 2 estimated region: {} - {}'.format(
			trans_ranges[1].lower, trans_ranges[1].upper))
	elif n_ranges == 1:
		print('Only one transition detected. Assuming start is replication site 1')
		print('Replication site 2 estimated region: {} - {}'.format(
			trans_ranges[0].lower, trans_ranges[0].upper))
	else:
		assert n_ranges == 0
		print('No inversions detected')

def load_data(out_dir):
	data_file = os.path.join(out_dir, 'matrix_data.npy')
	return np.load(data_file, allow_pickle=True)


def detect(out_dir, stride, window_len, padding, **args):
	ori_mat, window_locs, seq_len = load_data(out_dir)
	trans_idx = get_trans_idx(ori_mat)
	trans_ranges = get_trans_ranges(trans_idx, window_locs, window_len, stride, padding)
	print_transitions(trans_ranges, seq_len)
	return trans_ranges
