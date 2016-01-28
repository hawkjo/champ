import os
from collections import defaultdict


def phix_read_names(project_name, file_structure):
    fpath = os.path.join(
            file_structure.data_directory,
            'from_fourierseq',
            project_name,
            'phiX_mappings',
            'phiX_read_names.txt')
    return get_read_names(fpath)


def all_read_names(project_name, file_structure):
    fpath = os.path.join(
            file_structure.data_directory,
            'from_fourierseq',
            project_name,
            'all_fastqs',
            'all_read_names.txt')
    if not os.path.isdir(fpath):
        fpath = fpath.replace('all_fastqs', 'read_names')
    return get_read_names(fpath)


def get_read_names(fpath):
    with open(fpath) as f:
        tiles = defaultdict(set)
        for line in f:
            lane, tile = line.strip().rsplit(':', 4)[1:3]
            key = 'lane{0}tile{1}'.format(lane, tile)
            tiles[key].add(line.strip())
        return {key: list(values) for key, values in tiles.items()}
