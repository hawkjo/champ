"""
Variables and methods for finding local files and folders.
"""

import os
import numpy as np
from collections import defaultdict

base_dir = '/home/hawkjo/IlyaProjects/miseq_alignment'
fig_dir = os.path.join(base_dir, 'figs')
data_dir = os.path.join(base_dir, 'data')
fourier_data_dir = os.path.join(data_dir, 'from_fourierseq')
results_dir = os.path.join(base_dir, 'results')

targets = {
    'A': 'AAGGCCGAATTCTCACCGGCCCCAAGGTATTCAAG',
    'A-Csy': 'ACCGCCGAATTCTCACCGGCCCCAAGGTATTCAAG',
    'B': 'AAGTCGGCTCCTGTTTAGTTACGAGCGACATTGCT',
    'C': 'AAGCCAGTGATAAGTGGAATGCCATGTGGGCTGTC',
    'D': 'TTTAGTGATAAGTGGAATGCCATGTGG',
    'E': 'TTTAGACGCATAAAGATGAGACGCTGG'
}

def phiX_read_names_given_project_name(project_name):
    fpath = os.path.join(
            data_dir,
            'from_fourierseq',
            project_name,
            'phiX_mappings',
            'phiX_read_names.txt')
    return fastq_tiles_given_read_name_fpath(fpath)


def all_read_names_given_project_name(project_name):
    fpath = os.path.join(
            data_dir,
            'from_fourierseq',
            project_name,
            'all_fastqs',
            'all_read_names.txt')
    if not os.path.isdir(fpath):
        fpath = fpath.replace('all_fastqs', 'read_names')
    return fastq_tiles_given_read_name_fpath(fpath)


def fastq_tiles_given_read_name_fpath(fpath):
    fastq_tile_builder = defaultdict(set)
    for line in open(fpath):
        _, lane, tile, _, _ = line.strip().rsplit(':', 4)
        key = 'lane{0}tile{1}'.format(lane, tile)
        fastq_tile_builder[key].add(line.strip())
    return fastq_tiles_given_fastq_tile_builder(fastq_tile_builder)


def fastq_tiles_given_read_names(read_names):
    fastq_tile_builder = defaultdict(set)
    for read_name in read_names:
        _, lane, tile, _, _ = read_name.rsplit(':', 4)
        key = 'lane{0}tile{1}'.format(lane, tile)
        fastq_tile_builder[key].add(read_name)
    return fastq_tiles_given_fastq_tile_builder(fastq_tile_builder)

def fastq_tiles_given_fastq_tile_builder(fastq_tile_builder):
    fastq_tiles = {}
    for key, values in fastq_tile_builder.items():
        fastq_tiles[key] = list(values)
    return fastq_tiles
