#!/usr/bin/env python
import sys
import os
import argparse
from detect import orientation_matrix, eval_matrix, find_repeats, correct_inversion, uneven_regions
from detect.utils import utils, genome_utils


def check_args(args):
    chrom_len = genome_utils.chrom_len(args['genome_path'], chrom_id=args['chrom_id'])
    assert chrom_len/args['window_len'] > 2 and chrom_len/args['stride'] >2, \
            'window length and stride must be sufficiently large compared to chromosome length'

def add_args(args):
    out_path = lambda x: os.path.join(args['out_dir'], x)
    args['tmp_dir'] = out_path('tmp')
    tmp_path = lambda x: os.path.join(args['tmp_dir'], x)
    if args['chrom_id'] is None: 
        args['chrom_id'] = genome_utils.get_default_chrom_id(args['genome_path'])
    args['chrom_len'] = genome_utils.chrom_len(args['genome_path'])
    args['tmp_genome_path'] = tmp_path('genome.fasta')
    args['matrix_path'] = tmp_path('matrix_data.npy') 
    args['png_path'] = out_path('orientation_mat.png')
    args['repeats_path'] = tmp_path('nuc.coords')
    args['delta_path'] = tmp_path('nuc.delta')
    args['candidates_path'] = tmp_path('candidates.npy')
    args['corr_genome_path'] = utils.get_corrected_genome_path(args['genome_path'], args['out_dir'])
    args['corr_png_path'] = out_path('corrected_orientation_mat.png')
    args['corr_matrix_path'] = tmp_path('corr_matrix_data.npy')
    args['trans_data_path'] = tmp_path('trans_ranges.npy')
    args['summary'] = out_path('summary.txt')
    return args

def parse_args():
    parser = argparse.ArgumentParser(description='Detect inversions in a genome and correct misassemblies')
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
    parser.add_argument('--thresh', dest='len_thresh', type=int, default=10000,
            help='Minimum length threshold for candidate repeats')
    parser.add_argument('--genome', type=str, dest='genome_path', required=True,
            help='Path to input genome file. File types supported: .fasta, .fasta.gz')
    parser.add_argument('--out', type=str, dest='out_dir', required=True,
            help='Output directory')
    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parse_args()
    args = add_args(args)
    check_args(args)
    utils.setup_out_dir(**args)
    
    utils.print_banner('PARAMETERS')
    utils.print_matrix_info(**args)
    
    utils.print_banner('MAKING ORIENTATION MATRIX')
    orientation_matrix.make_matrix(**args)
    
    utils.print_banner('EVALUATING MATRIX')
    is_imbalanced = eval_matrix.evaluate(**args)

    if is_imbalanced:
        utils.print_banner('SEARCHING FOR REPEATS')
        print('RUNNING NUCMER')
        find_repeats.run_nucmer(**args)
        all_repeats = find_repeats.find_all_repeats(**args)
        
        candidates_found = uneven_regions.find_candidate_repeats(all_repeats, **args)
        if candidates_found:
            utils.print_banner('CORRECTING MISASSEMBLY')
            correct_inversion.correct(**args)
            
            utils.print_banner('MAKING CORRECTED MATRIX')
            orientation_matrix.make_corrected_matrix(**args)
            
            utils.print_banner('EVALUATING CORRECTED MATRIX')
            eval_matrix.evaluate_corrected(**args)
    
    if not args['save_tmp_files']:
        os.system('rm -r {}'.format(args['tmp_dir']))
        