import editdistance
import numpy as np
import random
import sys
from champ import initialize
import yaml


def get_max_edit_dist(target):
    dists = [editdistance.eval(target, rand_seq(target)) for _ in xrange(1000)]
    return min(10, np.percentile(dists, 0.5))


def rand_seq(target):
    seq_len = int(random.normalvariate(len(target), len(target) / 10))
    return ''.join(random.choice('ACGT') for _ in xrange(seq_len))


def get_target_reads(target, reads_by_seq_fpath, out_fpath):
    max_edit_dist = get_max_edit_dist(target)
    print('Max edit distance:', max_edit_dist)
    with open(out_fpath, 'w') as out:
        for line in open(reads_by_seq_fpath):
            words = line.strip().split()
            seq = words[0]
            read_names = words[1:]
            if editdistance.eval(target, seq) <= max_edit_dist:
                out.write('\n'.join(read_names) + '\n')


if __name__ == '__main__':
    usg_fmt = '{} <image_directory> <reads_by_seq_fpath> <out_fpath>'.format(sys.argv[0])
    if len(sys.argv) != len(usg_fmt.split()):
        sys.exit(usg_fmt)

    image_directory, reads_by_seq_fpath, out_fpath = sys.argv[1:]
    metadata = initialize.load(image_directory)
    with open("/shared/targets.yml") as f:
        target_sequences = yaml.load(f)

    get_target_reads(target_sequences[metadata['perfect_target_name']],
                     reads_by_seq_fpath,
                     out_fpath)
