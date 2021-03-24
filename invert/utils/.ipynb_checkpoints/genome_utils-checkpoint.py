import pandas as pd
import numpy as np
import numpy.random as rnd
import os
import gzip
from io import StringIO


DATA_PATH = '/home/a-m/gcgreen2/data/'
GENOMES_PATH = os.path.join(DATA_PATH,'genomic/')
BASES = {'A','T','C','G'}

assemblies_index = pd.read_csv(os.path.join(GENOMES_PATH, 'assemblies_index_full.txt'), sep='\t')
assemblies_index = assemblies_index.set_index('taxid')

def get_stringio(file):
    with gzip.open(file) as gfh: cov = gfh.read().decode('utf-8')
    return StringIO(cov)

def parse_fasta(genome):
    chromosomes = {}
    genome = genome.split('>')[1:]
    for chrom in genome:
        chrom = chrom.split('\n')
        chrom_id = chrom[0]
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
    chromosomes = parse_fasta(genome)
    return chromosomes

def get_chrom_number(chroms, number=0):
    return list(chroms.values())[number]

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
