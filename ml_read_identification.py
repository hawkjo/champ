import sys
import os
import random
import itertools
from collections import defaultdict
import numpy as np
from Bio import SeqIO
import pickle
import fastqtools
from adapters_cython import simple_hamming_distance
from misctools import gzip_friendly_open


bases = 'ACGT'
def rand_seq(seq_len):
    return ''.join(random.choice(bases) for _ in xrange(seq_len))


def get_max_ham_dists(min_len, max_len):
    dists = defaultdict(list)
    for _ in xrange(50000):
        ref_seq = rand_seq(max_len)
        new_seq = rand_seq(max_len)
        for i in range(min_len, max_len+1):
            dists[i].append(simple_hamming_distance(ref_seq[:i], new_seq[:i]))
    max_ham_dists = [min(np.percentile(dists[i], 0.1), int(i/4)) for i in range(min_len, max_len+1)]
    return max_ham_dists


def classify_reads(fastq_fpaths, log_p_fpath, min_len, max_len, out_fpath):
    """
    Classifies reads by overlapping ML sequence identity.
    """
    #--------------------------------------------------------------------------------
    # Load log_p dict of dicts of lists. Addessed as follows:
    #   
    #   log_p_struct[true_base][read_base][phred_score]
    #--------------------------------------------------------------------------------
    with open(log_p_fpath) as f:
        log_p_struct = pickle.load(f)

    #--------------------------------------------------------------------------------
    # Find max hamming distances per length considered
    #--------------------------------------------------------------------------------
    print 'Finding max hamming distances by length...',
    sys.stdout.flush()
    max_ham_dists = get_max_ham_dists(min_len, max_len)
    print 'Done'

    #--------------------------------------------------------------------------------
    # Make classifying function with given params
    #--------------------------------------------------------------------------------
    bases = 'ACGT'
    bases_set = set(bases)
    def classify_seq(rec1, rec2):
        # Store as strings
        seq1 = str(rec1.seq)
        seq2_rc = str(rec2.seq.reverse_complement())
        loc_max_len = min(max_len, len(seq1), len(seq2_rc))

        # Find aligning sequence, indels are not allowed, starts of reads included
        sig_lens = [i for i, max_ham in zip(range(min_len, loc_max_len + 1), max_ham_dists)
                    if simple_hamming_distance(seq1[:i], seq2_rc[-i:]) < max_ham]
        if len(sig_lens) != 1:
            return None

        seq2_len = sig_lens[0]
        seq2_match = seq2_rc[-seq2_len:]
        seq1_match = seq1[:seq2_len]

        # Get corresponding quality scores
        quals1 = rec1.letter_annotations['phred_quality'][:seq2_len]
        quals2 = rec2.letter_annotations['phred_quality'][::-1][-seq2_len:]

        # Build concensus sequence
        ML_bases = []
        for r1, q1, r2, q2 in zip(seq1_match, quals1, seq2_match, quals2):
            if r1 in bases and r1 == r2:
                ML_bases.append(r1)
            elif set([r1, r2]) <= bases_set and q1 > 2 and q2 > 2:
                r1_score = log_p_struct[r1][r1][q1] + log_p_struct[r1][r2][q2]
                r2_score = log_p_struct[r2][r1][q1] + log_p_struct[r2][r2][q2]
                if r1_score > r2_score:
                    ML_bases.append(r1)
                else:
                    ML_bases.append(r2)
            elif r1 in bases and q1 > 2:
                ML_bases.append(r1)
            elif r2 in bases and q2 > 2:
                ML_bases.append(r2)
            else:
                return None
        return ''.join(ML_bases)

    #--------------------------------------------------------------------------------
    # Pair fpaths and classify seqs
    #--------------------------------------------------------------------------------
    pe_fpaths, se_fpaths = fastqtools.find_paired_and_unpaired_files_from_fpaths(fastq_fpaths)
    assert not se_fpaths, 'All fastq files should be paired.'

    read_names_given_seq = defaultdict(list)
    for fpath1, fpath2 in pe_fpaths:
        print '{}, {}'.format(*map(os.path.basename, (fpath1, fpath2)))
        discarded = 0
        total = 0
        for i, (rec1, rec2) in enumerate(
                itertools.izip(SeqIO.parse(gzip_friendly_open(fpath1), 'fastq'),
                               SeqIO.parse(gzip_friendly_open(fpath2), 'fastq'))
        ):
            if i % 100000 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
            total += 1
            seq = classify_seq(rec1, rec2)
            if seq:
                read_names_given_seq[seq].append(str(rec1.id))
            else:
                discarded += 1
        found = total - discarded
        print
        print 'Found {} of {} ({:.1f}%)'.format(found, total, 100 * found / float(total))

    #--------------------------------------------------------------------------------
    # Output results
    #--------------------------------------------------------------------------------
    with open(out_fpath, 'w') as out:
        for seq, read_names in sorted(read_names_given_seq.items()):
            out.write('{}\t{}\n'.format(seq, '\t'.join(read_names)))


def isint(a):
    try:
        int(a)
        return float(a) == int(a)
    except:
        return False

if __name__ == '__main__':
    usage_fmt = '{} <min_len> <max_len> <max_mismatch> <out_fpath> <log_p_fpath> <fastq_fpaths>'.format(sys.argv[0])
    if len(sys.argv) < len(usage_fmt.split()):
        helpstr = """
Usage: {}

    min_len:            Minimum allowed overlap
    max_len:            Maximum allowed overlap
    out_fpath:          Location to write output file
    log_p_fpath:        Location of pickle file with probability struct
    fastq_fpaths:       List of all fastq files in run
""".format(usage_fmt)
        sys.exit(helpstr)

    min_len, max_len = map(int, sys.argv[1:3])
    assert all(map(isint, (min_len, max_len))), 'Min and Max lens and hams must be integers'
    out_fpath = sys.argv[3]
    log_p_fpath = sys.argv[4]
    fastq_fpaths = sys.argv[5:]

    classify_reads(fastq_fpaths, log_p_fpath, min_len, max_len, out_fpath)
