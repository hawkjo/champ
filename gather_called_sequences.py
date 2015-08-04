import os
import glob
import sys
from collections import defaultdict
from Bio import SeqIO
from misctools import gzip_friendly_open


def gather_called_sequences(fastq_dir, read_call_dir):
    print 'Reading called names'
    called_seq_names = set()
    for fpath in glob.glob(os.path.join(read_call_dir, '*_read_calls.txt')):
        called_seq_names.update(line.strip().split()[0] for line in open(fpath))
    print len(called_seq_names), 'called names'

    print 'Reading fastqs'
    out_seqs = defaultdict(list)
    for fpath in (glob.glob(os.path.join(fastq_dir, '*.fastq'))
                  + glob.glob(os.path.join(fastq_dir, '*.fq'))
                  + glob.glob(os.path.join(fastq_dir, '*.fastq.gz'))
                  + glob.glob(os.path.join(fastq_dir, '*.fq.gz'))):
        if '_I1_' in fpath or '_I2_' in fpath:
            continue
        for rec in SeqIO.parse(gzip_friendly_open(fpath), 'fastq'):
            if rec.id in called_seq_names:
                out_seqs[rec.id].append(rec)

    all_lengths = set(len(seq_list) for seq_list in out_seqs.values())
    assert len(all_lengths) == 1, all_lengths

    print 'Writing called sequences'
    out_fpath = os.path.join(read_call_dir, 'gathered_sequences.fastq')
    with open(out_fpath, 'w') as out:
        for name, recs in out_seqs.items():
            SeqIO.write(recs, out, 'fastq')


if __name__ == '__main__':
    arg_fmt = '{0} <fastq_dir> <read_call_dir>'.format(sys.argv[0])
    if len(sys.argv) != len(arg_fmt.split()):
        sys.exit('Usage: {0}'.format(arg_fmt))

    fastq_dir = sys.argv[1]
    read_call_dir = sys.argv[2]

    gather_called_sequences(fastq_dir, read_call_dir)
