import sys
import os
import argparse
from invert import orientation_matrix, detect_inversions, find_repeats, correct_inversion
from invert.utils import genome_utils

def print_banner(x):
	print('\n=================== {} ==================='.format(str(x)))

def print_matrix_info(genome_file, window_len, stride, out_dir, chrom_id, **args):
	print('genome file: {}'.format(genome_file))
	print('sequence id chosen: {}'.format(chrom_id))
	print('window length: {},  stride: {}'.format(window_len, stride))
	print('output directory: {}'.format(out_dir))

def setup_out_dir(out_dir, genome_file, **args):
	os.makedirs(os.path.join(out_dir, 'tmp'), exist_ok=True)
	tmp_genome_file = os.path.join(out_dir,'tmp','genome.fasta')
	chroms = genome_utils.get_chromosomes(genome_file)
	genome_utils.write_fasta(tmp_genome_file, chroms)
	return out_dir, tmp_genome_file

def get_default_chrom_id(genome_file):
	chroms = genome_utils.get_chromosomes(genome_file)
	return list(chroms.keys())[0]

def get_corrected_args(cor_genome_file, **args):
	corrected_args = args.copy()
	corrected_args['save_png'] = True
	corrected_args['png_name'] = 'orientation_mat_corrected.png'
	corrected_args['genome_file'] = cor_genome_file
	return corrected_args

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
	parser.add_argument('--save-files', action='store_true', dest='save_tmp_files', default=False,
			help='Save temporary files (nucmer outputs, unzipped genome, matrix data)')
	parser.add_argument('--nucmer', dest='nucmer_path', type=str, default='nucmer',
			help='Path to nucmer executable (default assumes already in PATH)')
	parser.add_argument('--genome', type=str, dest='genome_file', required=True,
			help='Genome file. File types supported: .fasta, .fasta.gz')
	parser.add_argument('--out', type=str, dest='out_dir', required=True,
			help='Output directory')
	args = vars(parser.parse_args())
	args['png_name'] = 'orientation_mat.png'
	return args

if __name__ == "__main__":
	args = parse_arguments()
	if args['chrom_id'] is None: args['chrom_id'] = get_default_chrom_id(args['genome_file'])
	out_dir, tmp_genome_file = setup_out_dir(**args)
	
	print_banner('PARAMETERS')
	print_matrix_info(**args)
	
	print_banner('MAKING ORIENTATION MATRIX')
	orientation_matrix.make_matrix(**args)

	print_banner('DETECTING INVERSIONS')
	trans_ranges = detect_inversions.detect(**args)
	if args['correct_inv']:
		print_banner('RUNNING NUCMER')
		os.system('{0} -r -l 15 -c 200 -g 20 -o -p {1}/tmp/nuc {2} {2}'.format(
			args['nucmer_path'], out_dir, tmp_genome_file))
		
		print_banner('SEARCHING FOR REPEATS')
		best_repeat = find_repeats.find_best(trans_ranges=trans_ranges, **args)

		if best_repeat is not None:
			print_banner('CORRECTING MISASSEMBLY')
			cor_genome_file = correct_inversion.correct(best_repeat=best_repeat, **args)
			print_banner('MAKING CORRECTED ORIENTATION MATRIX')
			corrected_args = get_corrected_args(cor_genome_file, **args)
			orientation_matrix.make_matrix(**corrected_args)
	
	if not args['save_tmp_files']:
		os.system('rm -r {}'.format(os.path.join(out_dir,'tmp')))
		
