import editdistance
import numpy as np
import random
import sys
import os

targets = {
    'A': 'AAGGCCGAATTCTCACCGGCCCCAAGGTATTCAAG',
    'A-Csy': 'ACCGCCGAATTCTCACCGGCCCCAAGGTATTCAAG',
    'B': 'AAGTCGGCTCCTGTTTAGTTACGAGCGACATTGCT',
    'C': 'AAGCCAGTGATAAGTGGAATGCCATGTGGGCTGTC',
    'D': 'TTTAGTGATAAGTGGAATGCCATGTGG',
    'E': 'TTTAGACGCATAAAGATGAGACGCTGG'
}


def get_max_edit_dist(target):
    dists = [editdistance.eval(target, rand_seq(target)) for _ in xrange(1000)]
    return min(10, np.percentile(dists, 0.5))


def rand_seq(target):
    seq_len = int(random.normalvariate(len(target), len(target) / 10))
    return ''.join(random.choice('ACGT') for _ in xrange(seq_len))


def get_target_reads(target, reads_by_seq_fpath):
    max_edit_distance = get_max_edit_dist(target)
    print('Max edit distance:', max_edit_distance)
    for line in open(reads_by_seq_fpath):
        words = line.strip().split()
        seq = words[0]
        read_names = words[1:]
        distance = editdistance.eval(target, seq)
        if distance <= max_edit_distance:
            yield '\n'.join(read_names) + '\n'


if __name__ == '__main__':
    usg_fmt = '{} <reads_by_seq_fpath> <out_dir>'.format(sys.argv[0])
    if len(sys.argv) != len(usg_fmt.split()):
        sys.exit(usg_fmt)

    reads_by_seq_fpath, out_dir = sys.argv[1:]
    for label, target_sequence in targets.items():
        out_path = os.path.join(out_dir, "target_{}_read_names.txt".format(label.lower()))
        with open(out_path, 'w+') as out:
            for reads in get_target_reads(target_sequence, reads_by_seq_fpath):
                out.write(reads)
