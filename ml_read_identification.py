import sys
import os
import itertools
from collections import defaultdict
from Bio import SeqIO
import pickle
import fastqtools
from adapters_cython import simple_hamming_distance
from misctools import gzip_friendly_open


def classify_reads(fastq_fpaths, log_p_fpath, min_len, max_len, max_ham, out_fpath):
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
    # Make classifying function with given params
    #--------------------------------------------------------------------------------
    bases = 'ACGT'
    bases_set = set(bases)
    def classify_seq(rec1, rec2):
        # Store as strings
        seq1 = str(rec1.seq)
        seq2_rc = str(rec2.seq.reverse_complement())

        # Find aligning sequence, indels are not allowed, starts of reads included
        hams = [simple_hamming_distance(seq1[:i], seq2_rc[-i:])
                for i in range(min_len, max_len + 1)]
        min_ham = min(hams)
        if min_ham > max_ham:
            return None

        seq2_len = min(range(min_len, max_len + 1), key=lambda i: hams[i-min_len])
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
        print
        print 'Discarded {} of {} ({:.1f}%)'.format(discarded, total, 100 * discarded / float(total))

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
    max_mismatch:       Maximum allowed mismatch between reads
    out_fpath:          Location to write output file
    log_p_fpath:        Location of pickle file with probability struct
    fastq_fpaths:       List of all fastq files in run
""".format(usage_fmt)
        sys.exit(helpstr)

    assert all(map(isint, sys.argv[1:4])), 'Min and Max lens and hams must be integers'
    min_len, max_len, max_ham = map(int, sys.argv[1:4])
    out_fpath = sys.argv[4]
    log_p_fpath = sys.argv[5]
    fastq_fpaths = sys.argv[6:]

    classify_reads(fastq_fpaths, log_p_fpath, min_len, max_len, max_ham, out_fpath)
