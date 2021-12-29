import sys
import os
import numpy as np

from detect.utils import genome_utils

def print_banner(x):
	print('\n=================== {} ==================='.format(str(x)))

def print_header():
	print_banner('ORIENTED - TNF')
	for _ in range(3): print('=======================================')
		
def print_matrix_info(genome_path, window_len, stride, out_dir, chrom_id, summary, **args):
	info = 'genome file: {}\n'.format(genome_path) + \
			'sequence id chosen: {}\n'.format(chrom_id) + \
			'window length: {},  stride: {}\n'.format(window_len, stride) + \
			'output directory: {}\n\n'.format(out_dir)
	print(info)
	with open(summary, 'a') as fh:
		fh.write(info)
	
def setup_out_dir(out_dir, genome_path, tmp_genome_path, **args):
	tmp_dir = os.path.join(out_dir, 'tmp')
	os.makedirs(tmp_dir, exist_ok=True)
	chroms = genome_utils.get_chromosomes(genome_path)
	genome_utils.write_fasta(tmp_genome_path, chroms)
	
def get_corrected_genome_path(genome_path, out_dir, **args):
# 	index = '' if index==0 else str(index+1)
	corrected_path = os.path.split(genome_path)[1]
	prefix, suffix = os.path.splitext(corrected_path)
	if suffix == '.gz':
		prefix, suffix = os.path.splitext(prefix)
	if suffix == '': suffix = 'fasta'
	corrected_path = prefix + '_inv' + suffix
	return os.path.join(out_dir, corrected_path)

def get_corrected_png_path(out_dir, **args):
#	index = '' if index==0 else str(index+1)
	return os.path.join(out_dir, 'corrected_orientation_mat.png')

def save_file(filepath, x):
	x = np.array(x, dtype=object)
	np.save(filepath, x, allow_pickle=True)
	
def load_file(filepath):
	return np.load(filepath, allow_pickle=True)
