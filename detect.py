import sys
import os
import argparse
from invert import orientation_matrix
from invert import detect_inversions
from invert.utils import genome_utils

def print_banner(x):
	print('=================== {} ==================='.format(str(x)))

def print_matrix_info(genome_file,window_len,stride,out_dir,chrom_id,**kwargs):
	chroms = genome_utils.get_chromosomes(genome_file)
	if chrom_id is None: chrom_id = list(chroms.keys())[0]
	print('genome file: {}'.format(genome_file))
	print('sequence id chosen: {}'.format(chrom_id))
	print('window length: {},  stride: {}'.format(window_len, stride))
	print('output directory: {}'.format(out_dir))

def parse_arguments():
	parser = argparse.ArgumentParser(description='Detect inversions in a genome')
	parser.add_argument('--chrom', type=str, dest='chrom_id',
			help='Chromsome/Sequence id')
	parser.add_argument('--window', type=int, default=50000, dest='window_len', 
			help='Window length for matrix calculation. Default: 5e4')
	parser.add_argument('--stride', type=int, default=25000,
			help='Stride between windows. Default: 2.5e4')
	parser.add_argument('--pad', type=int, default=0, dest='padding',
			help='Padding for repeat search range. Default: 0')
	parser.add_argument('--no-correct', action='store_false', dest='correct_inv', default=True, 
			help='Do not correct detected misassembly')
	parser.add_argument('--no-png', action='store_false', dest='save_png', default=True,
			help='Do not save orientation matrix heatmap png.')
	parser.add_argument('genome_file', type=str,
			help='Genome file. File types supported: .fasta, .fasta.gz, .fastq, .fastq.gz')
	parser.add_argument('out_dir', type=str,
			help='Output directory')
	args = vars(parser.parse_args())
	return args

if __name__ == "__main__":
	args = parse_arguments()
	print_matrix_info(**args)
	os.makedirs(args['out_dir'], exist_ok=True)
	orientation_matrix.make_matrix(**args)
	ranges = detect_inversions.detect(**args)
	#if correct_inv:
		

