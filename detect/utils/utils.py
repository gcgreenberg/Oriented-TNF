import sys
import os
import numpy as np

from detect.utils import genome_utils

def print_banner(x):
	print('\n=================== {} ==================='.format(str(x)))

def print_header():
	print_banner('ORIENTED - TNF')
	for _ in range(3): print('=======================================')
		
def print_matrix_info(genome_path, window_len, stride, out_dir, chrom_id, **args):
	print('genome file: {}'.format(genome_path))
	print('sequence id chosen: {}'.format(chrom_id))
	print('window length: {},  stride: {}'.format(window_len, stride))
	print('output directory: {}'.format(out_dir))
	
def setup_out_dir(out_dir, genome_path, tmp_genome_path, **args):
	tmp_dir = os.path.join(out_dir, 'tmp')
	os.makedirs(tmp_dir, exist_ok=True)
	chroms = genome_utils.get_chromosomes(genome_path)
	genome_utils.write_fasta(tmp_genome_path, chroms)
	
def get_corrected_genome_path(genome_path, out_dir, index=0, **args):
	index = '' if index==0 else str(index+1)
	corrected_path = os.path.split(genome_path)[1]
	prefix, suffix = os.path.splitext(corrected_path)
	if suffix == '.gz':
		prefix, suffix = os.path.splitext(prefix)
	if suffix == '': suffix = 'fasta'
	corrected_path = prefix + '_inv' + index + suffix
	return os.path.join(out_dir, corrected_path)

def get_corrected_png_path(out_dir, index=0, **args):
	index = '' if index==0 else str(index+1)
	return os.path.join(out_dir, 'corrected_orientation_mat{}.png'.format(index))

def save_file(filepath, x):
	x = np.array(x, dtype=object)
	np.save(filepath, x, allow_pickle=True)
	
def load_file(filepath):
	return np.load(filepath, allow_pickle=True)