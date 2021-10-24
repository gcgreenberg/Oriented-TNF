import os
import gzip
import pandas as pd
import numpy as np

BASES = {'A','T','C','G'}

def parse_fasta(genome):
	chromosomes = {}
	genome = genome.split('>')[1:]
	for chrom in genome:
		chrom = chrom.split('\n')
		chrom_id = chrom[0].split(' ')[0]
		chrom = ''.join(chrom[1:])
		for char in set(chrom).difference({'A','T','C','G'}):
			chrom = chrom.replace(char,'')
		chromosomes[chrom_id] = chrom
	return chromosomes
	
def get_chromosomes(file):
	if os.path.splitext(file)[1] == '.gz':
		with gzip.open(file) as fh:
			genome = fh.read().decode('utf-8').rstrip()
	else:
		with open(file) as fh:
			genome = fh.read().rstrip()
	chroms = parse_fasta(genome)
	return chroms

def get_chromosome(file, chrom_id=None, largest=False):
	chroms = get_chromosomes(file)
	if chrom_id is None:
		return list(chroms.values())[0]
	else:
		return chroms[chrom_id]
	
def get_default_chrom_id(genome_file):
	chroms = get_chromosomes(genome_file)
	return list(chroms.keys())[0]

def get_concat_genome(file):
	chroms = get_chroms(file)
	genome = ''.join(chroms.values())
	return genome

########### WRITING FILES

def write_fasta(file, chroms):
	with open(file, 'w') as fh:
		for seqid in chroms:
			fh.write('>' + seqid + "\n")
			fh.write(chroms[seqid] + "\n")
