"""
We wish to find all phiX reads to which a given oligo will bind. This requires strand-sensitive
inference of overlap with the sequence in a given read.
"""

import glob
import os
import sys
from collections import Counter

import pysam
from Bio import SeqIO
from Bio.Seq import Seq


def find_oligo_sites(oligo_seq, phiX_genome_fpath):
    recs = list(SeqIO.parse(phiX_genome_fpath, 'fasta'))
    assert len(recs) == 1, 'PhiX genome must be one sequence.'
    phiX_seq = str(recs[0].seq)
    phiX_seq_rc = str(recs[0].seq.reverse_complement())
    fwd_sites, rc_sites = [], []
    # Forward sites
    i = 0
    while True:
        try:
            site = i + phiX_seq[i:].index(oligo_seq)
            fwd_sites.append(site)
            i = site + 1
        except ValueError:
            break
    # RC sites
    i = 0
    while True:
        try:
            site = i + phiX_seq_rc[i:].index(oligo_seq)
            rc_site_in_fwd = len(phiX_seq) - site - len(oligo_seq)
            assert str(Seq(phiX_seq[rc_site_in_fwd:rc_site_in_fwd + len(oligo_seq)]).reverse_complement()) == oligo_seq
            rc_sites.append(rc_site_in_fwd)
            i = site + 1
        except ValueError:
            break
    return fwd_sites, rc_sites


def infer_reads(oligo_seq, phiX_genome_fpath, phiX_mapping_dir, out_fpath):
    """
    Infer which phiX reads the given oligo sequence will bind to.
    """
    # To infer reads to which an oligo will bind, we need the following picture of ssDNA strand
    # attached to an illumina chip after completion of sequencing:
    #
    #
    #               3'  |
    #            R2_rc  |
    #                   |  *
    #                   | /
    #                   ||   Oligo
    #                   ||
    #               R1  |
    #               5'  |
    # ---------------------------------------
    #
    # So in mapped space:
    #
    #       PhiX genome: --- R1> ---- <R2 ----
    #
    # means we should be looking for the oligo in the reverse complement of phiX and
    #
    #       PhiX genome: --- R2> ---- <R1 ----
    #
    # means we should be looking for the oligo in the forwards sequence.

    oligo_seq_rc = str(Seq(oligo_seq).reverse_complement())
    bam_fpaths = glob.glob(os.path.join(phiX_mapping_dir, '*.bam'))
    out_read_names = set()  # Use set structure due to multiple mappings
    properly_paired_read_names = set()
    unpaired_read_names = set()
    stats = Counter()
    fwd_sites, rc_sites = find_oligo_sites(oligo_seq, phiX_genome_fpath)
    # fwd_sites, rc_sites = rc_sites, fwd_sites
    print '{} forward site(s), {} Reverse complement site(s)'.format(len(fwd_sites), len(rc_sites))

    for bam_fpath in bam_fpaths:
        for read in pysam.Samfile(bam_fpath):
            if read.is_proper_pair:
                # Determine orientation and appropriate strand for oligo
                if read.isize > 0:
                    # read is the "left" read
                    start = read.pos
                    end = start + read.isize
                    if read.is_read1:
                        sites = rc_sites
                        stats['Properly paired left R1 read'] += 1
                    else:
                        sites = fwd_sites
                        stats['Properly paired left R2 read'] += 1
                else:
                    # read is the "right" read
                    start = read.mpos
                    end = start + abs(read.isize)
                    if read.is_read1:
                        sites = fwd_sites
                        stats['Properly paired right R1 read'] += 1
                    else:
                        sites = rc_sites
                        stats['Properly paired right R2 read'] += 1
                # Test for oligo
                for site in sites:
                    if start <= site and site + len(oligo_seq) <= end:
                        out_read_names.add(read.qname)
                        properly_paired_read_names.add(read.qname)
                        break  # Just for efficiency
            else:
                # Read is not properly paired. Just look in sequence.
                if read.is_read1:
                    stats['Unpaired R1 read'] += 1
                    if oligo_seq_rc in read.seq:
                        out_read_names.add(read.qname)
                        unpaired_read_names.add(read.qname)
                else:
                    stats['Unpaired R2 read'] += 1
                    if oligo_seq in read.seq:
                        out_read_names.add(read.qname)
                        unpaired_read_names.add(read.qname)

    print 'Properly paired seqs found: {}'.format(len(properly_paired_read_names))
    print 'Unpaired reads found: {}'.format(len(unpaired_read_names))
    for stat, count in sorted(stats.items()):
        print '{} reads total: {}'.format(stat, count)

    with open(out_fpath, 'w') as out:
        out.write('\n'.join(out_read_names))


if __name__ == '__main__':
    usg_fmt = '{} <oligo_seq> <phiX_genome_fasta> <phiX_mapping_dir> <out_fpath>'.format(sys.argv[0])
    if len(sys.argv) != len(usg_fmt.split()):
        sys.exit('Usage: ' + usg_fmt)

    infer_reads(*sys.argv[1:])
