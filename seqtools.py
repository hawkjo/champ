from collections import defaultdict
import sys
from champ.adapters_cython import simple_hamming_distance


def build_read_names_given_seq(target,
                               read_names_by_seq_fpath,
                               allowed_read_names_set,
                               is_interesting_seq,
                               max_ham,
                               verbose=True):
    interesting_reads = defaultdict(set)
    i = 0
    for i, line in enumerate(open(read_names_by_seq_fpath)):
        if verbose and i % 10000 == 0:
            sys.stdout.write('.')
            sys.stdout.flush()

        words = line.strip().split()
        seq = words[0]
        if is_interesting_seq(seq):
            read_names = set(words[1:]) & allowed_read_names_set
            interesting_reads[seq].update(read_names)
            last_start = len(seq) - len(target)
            if last_start < 0:
                continue
            min_ham_idx = min(range(0, last_start+1),
                              key=lambda i: simple_hamming_distance(target, seq[i:i+len(target)]))
            min_ham = simple_hamming_distance(target, seq[min_ham_idx:min_ham_idx+len(target)])
            if min_ham <= max_ham:
                min_ham_seq = seq[min_ham_idx:min_ham_idx+len(target)]
                interesting_reads[min_ham_seq].update(read_names)
    return interesting_reads