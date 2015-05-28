"""
Variables and methods for finding local files and folders.
"""

import os
import numpy as np
from collections import defaultdict

base_dir = '/home/hawkjo/BillProjects/constellation'
fig_dir = os.path.join(base_dir, 'figs')
data_dir = os.path.join(base_dir, 'data')
jah_base_dir = '/home/jah/projects/ilya/experiments'

def phiX_read_names_given_project_name(pname):
    fpath = os.path.join(
            data_dir,
            'from_fourierseq',
            pname,
            'phiX_mappings',
            'phiX_read_names.txt')

    fastq_tile_builder = defaultdict(set)
    for line in open(fpath):
        _, lane, tile, _, _ = line.strip().rsplit(':', 4)
        key = 'lane{0}tile{1}'.format(lane, tile)
        fastq_tile_builder[key].add(line.strip())

    fastq_tiles = {}
    for key, values in fastq_tile_builder.items():
        fastq_tiles[key] = list(values)
    return fastq_tiles
