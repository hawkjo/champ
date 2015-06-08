import os
import glob
import sys
import pysam
import itertools
from collections import defaultdict, Counter
import re
from Bio import SeqIO
from general_sequence_tools import dna_rev_comp
from adapters import build_adapters
from adapters_cython import simple_hamming_distance
from misctools import gzip_friendly_open


def get_restriction_site_info(r_site, bam_dir, out_fpath):
    r_site = r_site.upper()
    r_site_rc = dna_rev_comp(r_site)
    def r_site_locations(dna_seq):
        dna_seq = dna_seq.upper()
        return [m.start() for m in itertools.chain(re.finditer(r_site, dna_seq),
                                                   re.finditer(r_site_rc, dna_seq))]

    phix_double_seq = str(next(SeqIO.parse(open('/home/hawkjo/genomes/phix/phiX_double.fa'), 'fasta')).seq)

    adapter_in_R1, adapter_in_R2 = build_adapters(primer_type='PE')
    phiX_insert_len_est = 400

    r_site_positions = defaultdict(set)

    read_info_titles = ['qname', 'pos', 'qstart', 'qend', 'qlen', 'is_reverse', 'alen', 'tlen', 'seq']

    stats = Counter()
    bam_fpaths = glob.glob(os.path.join(bam_dir, '*.bam'))
    for bam_fpath in bam_fpaths:
        for read in pysam.Samfile(bam_fpath):
            if read.is_proper_pair and 100 < abs(read.tlen) < 1000:
                if read.tlen > 0 and read.is_reverse:
                    stats['proper pair, pos tlen/is_reverse=True'] += 1
                elif read.tlen <= 0 and not read.is_reverse:
                    stats['proper pair, neg tlen/is_reverse=False'] += 1
                elif read.is_read2:
                    continue
                elif read.is_read1:
                    # We recreate an approximation of the read, with the simplifying assumption
                    # that phiX goes all the way to the beginning of read2.
                    info_str = 'fpath=%s\n%s' \
                            % (bam_fpath, ' '.join('{0}={1}'.format(t, getattr(read, t)) for t in read_info_titles))
                    if read.tlen > 0:
                        stats['proper pair, tlen > 0'] += 1
                        read_seq = read.seq
                        template_seq = phix_double_seq[read.pos:read.pos + read.tlen]
                    else:
                        stats ['proper pair, tlen <= 0'] += 1
                        read_seq = dna_rev_comp(read.seq)
                        start = read.pos + read.qend + read.tlen
                        end = read.pos + read.qend
                        template_seq = dna_rev_comp(phix_double_seq[start:end])
                    assert len(template_seq) < 1000
                    if not simple_hamming_distance(read_seq[read.qstart:read.qend],
                                                   template_seq[:read.qlen]) < 0.3 * read.alen:
                        stats['bad hamming'] += 1
    
                    dna_seq = (dna_rev_comp(adapter_in_R2) 
                               + read.seq[:read.qstart]
                               + template_seq
                               + adapter_in_R1)  
                    r_site_positions[read.qname].update(r_site_locations(dna_seq))
                    continue
                else:
                    raise ValueError('read neither read1 or read2')
            else:
                stats['improper pair'] += 1

            # All reads that were not properly paired in a reconstructable manner end up here. We
            # add all positions based strictly on the read information, using the
            # phiX_insert_len_est to estimate positions of restriction sites in R2.
            if read.is_read1:
                dna_seq = dna_rev_comp(adapter_in_R2) + read.seq
                r_site_positions[read.qname].update(r_site_locations(dna_seq))
            else:
                dna_seq = dna_rev_comp(adapter_in_R1) + read.seq
                positions = [len(adapter_in_R2) + phiX_insert_len_est + len(adapter_in_R1) - pos
                            for pos in r_site_locations(dna_seq)]
                r_site_positions[read.qname].update(positions)

    for k, v in stats.items():
        print k, v
    with open(out_fpath, 'w') as out:
        info_iter = ((name, ','.join(str(pos) for pos in sorted(positions)))
                     for name, positions in r_site_positions.items())

        out.write('\n'.join('{0}\t{1}'.format(name, pos_str) for name, pos_str in info_iter))


if __name__ == '__main__':
    input_fmt = '{0} <restriction_site> <bam_dir> <out_file>'.format(sys.argv[0])
    if len(sys.argv) != len(input_fmt.split()):
        sys.exit('Usage: {0}'.format(input_fmt))
    get_restriction_site_info(*sys.argv[1:])
