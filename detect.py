import sys
import os
import argparse
from invert import orientation_matrix


def print_banner(x):
	print('=================== {} ==================='.format(str(x)))

def print_matrix_info(genome_file,window_len,stride,out_dir,chrom_id,**kwargs):
	#chroms = genome_utils.get_chromosomes(file)
	print('file: {}'.format(genome_file))
	#print('gen len: {}, num chrom: {}'.format(len(genome_utils.get_chrom_number(chroms,0)), len(chroms)))
	print('c_len: {},  stride: {}'.format(window_len, stride))
	print('out directory: {}'.format(out_dir))
	#seq_id = list(chroms.keys())[chrom_num]
	#seq_id = seq_id.split(' ')[0]
	print('seq_id: {}'.format(chrom_id))

def parse_arguments():
	parser = argparse.ArgumentParser(description='Detect inversions in a genome')
	parser.add_argument('--no-png', action='store_false', dest='save_png',
			help='Do not save orientation matrix heatmap png.')
	parser.add_argument('--chrom', type=str, dest='chrom_id',
			help='Chromsome/Sequence id')
	parser.add_argument('--window', type=int, default=50000, dest='window_len', 
			help='Window length for matrix calculation. Default: 5e4')
	parser.add_argument('--stride', type=int, default=25000,
			help='Stride between windows. Default: 2.5e4')
	parser.add_argument('--pad', type=int, default=0, dest='padding',
			help='Padding for repeat search range. Default: 0')
	parser.add_argument('genome_file', type=str,
			help='Genome file. File types supported: .fasta, .fasta.gz, .fastq, .fastq.gz')
	parser.add_argument('out_dir', type=str,
			help='Output directory')
	args = vars(parser.parse_args())
	return args

def test(**args):
	print(args)

if __name__ == "__main__":
	args = parse_arguments()
	print(args)
	print_matrix_info(**args)
	if not os.path.exists(args['out_dir']): 
		os.makedirs(args['out_dir'])
	orientation_matrix.make_matrix(**args)
