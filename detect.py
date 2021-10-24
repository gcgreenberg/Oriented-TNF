import sys
import os
import argparse
from detect import orientation_matrix, detect_inversions, find_repeats, correct_inversion
from detect.utils import utils, genome_utils


def add_args(args):
	out_dir = args['out_dir']
	args['tmp_dir'] = os.path.join(out_dir, 'tmp')
	tmp_path = lambda x: os.path.join(args['tmp_dir'], x)
	if args['chrom_id'] is None: 
		args['chrom_id'] = genome_utils.get_default_chrom_id(args['genome_path'])
	args['tmp_genome_path'] = tmp_path('genome.fasta')
	args['matrix_data_path'] = tmp_path('matrix_data.npy') 
	args['png_path'] = os.path.join(out_dir, 'orientation_mat.png')
	args['repeats_path'] = tmp_path('nuc.coords')
	args['delta_path'] = tmp_path('nuc.delta')
	args['candidates_path'] = tmp_path('candidates.npy')
	args['trans_data_path'] = tmp_path('trans_ranges.npy')
	return args

def parse_args():
	parser = argparse.ArgumentParser(description='Detect inversions in a genome')
	parser.add_argument('--chrom', type=str, dest='chrom_id',
			help='Chromsome/Sequence id')
	parser.add_argument('--window', type=int, default=100000, dest='window_len', 
			help='Window length for matrix calculation. Default: 1e5')
	parser.add_argument('--stride', type=int, default=25000,
			help='Stride between windows. Default: 2.5e4')
	parser.add_argument('--pad', type=int, default=0, dest='padding',
			help='Padding for repeat search range. Default: 0')
	parser.add_argument('--save-files', action='store_true', dest='save_tmp_files',
			help='Save temporary files (nucmer outputs, unzipped genome, matrix data)')
	parser.add_argument('--mummer', dest='mummer_path', type=str, default='',
			help='Path to mummer executables directory (default assumes already in PATH)')
	parser.add_argument('--genome', type=str, dest='genome_path', required=True,
			help='Path to input genome file. File types supported: .fasta, .fasta.gz')
	parser.add_argument('--out', type=str, dest='out_dir', required=True,
			help='Output directory')
	return vars(parser.parse_args())

if __name__ == "__main__":
	args = parse_args()
	args = add_args(args)
	utils.setup_out_dir(**args)
	
	utils.print_banner('PARAMETERS')
	utils.print_matrix_info(**args)
	
	utils.print_banner('MAKING ORIENTATION MATRIX')
	orientation_matrix.make_matrix(**args)
	
	utils.print_banner('DETECTING INVERSIONS')
	inversion_detected = detect_inversions.detect(**args)

	if inversion_detected:
		utils.print_banner('RUNNING NUCMER')
		find_repeats.run_nucmer(**args)
		
		utils.print_banner('SEARCHING FOR REPEATS')
		candidates_found = find_repeats.find_candidate_repeats(**args)

		if candidates_found:
			correct_inversion.correct(**args)
	
	if not args['save_tmp_files']:
		os.system('rm -r {}'.format(args['tmp_dir']))
		
