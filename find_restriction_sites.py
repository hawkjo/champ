import sys
import os
import glob
import re
from Bio import SeqIO
from collections import defaultdict, Counter
import local_config
from misctools import gzip_friendly_open
from adapters import build_adapters
from fastqtools import collapse_if_overlapped_pair, find_paired_and_unpaired_files
from general_sequence_tools import dna_rev_comp
import itertools


def make_restriction_site_position_finder(rest_site, exist_paired, assumed_gap_size=0):
    # Adapters
    ainr1, ainr2 = build_adapters()
    a_pre_r1 = dna_rev_comp(ainr2)
    a_pre_r2 = dna_rev_comp(ainr1)

    # Restriction site regex
    if rest_site == dna_rev_comp(rest_site):  # if palendromic
        pattern = rest_site
    else:
        pattern = '(%s|%s)' % (rest_site, dna_rev_comp(rest_site))
    rest_re = re.compile(pattern)

    def SE_restriction_sites(reads_dict):
        assert set(1) <= set(reads_dict.keys()), reads_dict.keys()
        if 1 in reads_dict:
            return [m.start() for m in rest_re.finditer(a_pre_r1 + reads_dict[1])]
        else:
            return []

    def PE_restriction_sites(reads_dict):
        assert set([1, 2]) <= set(reads_dict.keys()), reads_dict.keys()
        if 1 in reads_dict and 2 in reads_dict:  # is paired
            read_obj = collapse_if_overlapped_pair([reads_dict[1], reads_dict[2]])
            if isinstance(read_obj, str):
                return [m.start() for m in rest_re.finditer(a_pre_r1 + read_obj + ainr1)]
            else:
                assert len(read_obj) == 2, read_obj
                read_obj = [a_pre_r1 + read_obj[0], a_pre_r2 + read_obj[1]]
                read1_pos = [m.start() for m in rest_re.finditer(read_obj[0])] 
                read2_pos = [len(read_obj[0]) + assumed_gap_size + len(read_obj[1])  # total len
                             - m.start() - len(rest_site)  # minus that known to be to the right
                             for m in rest_re.finditer(read_obj[1])]
                return read1_pos + read2_pos
        elif 1 in reads_dict:
            return [m.start() for m in rest_re.finditer(a_pre_r1 + reads_dict[1])]
        elif 2 in reads_dict:
            read_obj = a_pre_r2 + reads_dict[2]
            return [2*len(read_obj) + assumed_gap_size - m.start() - len(rest_site)
                    for m in rest_re.finditer(read_obj[1])]
        else:
            return []

    if exist_paired:
        return PE_restriction_sites
    else:
        return SE_restriction_sites


def find_restriction_sites(project_name, rest_name, rest_site):
    fq_dir = os.path.join(local_config.fourier_data_dir, project_name, 'all_fastqs')
    # --------------------------------------------------------------------------------
    # Gather all read (pairs)
    # --------------------------------------------------------------------------------
    print 'Loading reads'
    read_pairs = defaultdict(dict)
    exist_paired = False
    for fpath in glob.glob(os.path.join(fq_dir, '*.fastq*')):
        fname = os.path.basename(fpath)
        if '_R1' in fname:
            read_idx = 1
        elif '_R2' in fname:
            read_idx = 2
            exist_paired = True
        else:
            continue

        print fname
        for rec in SeqIO.parse(gzip_friendly_open(fpath), 'fastq'):
            assert read_idx not in read_pairs[rec.id], rec.description
            read_pairs[rec.id][read_idx] = str(rec.seq)

    # --------------------------------------------------------------------------------
    # Find restriction sites
    # --------------------------------------------------------------------------------
    print('Finding restriction sites')
    restriction_site_positions = make_restriction_site_position_finder(rest_site, exist_paired)
    out_fpath = os.path.join(fq_dir, rest_name + '_restriction_sites.txt')
    with open(out_fpath, 'w') as out:
        for read_name, reads_dict in read_pairs.items():
            positions = restriction_site_positions(reads_dict)
            out.write('%s\t%s\n' % (read_name, ','.join('%d' % pos for pos in positions)))


def count_restriction_sites(project_name, rest_name, rest_site, max_reads_per_file=None):
    fq_dir = os.path.join(local_config.fourier_data_dir, project_name, 'all_fastqs')
    pe_fpaths, se_fpaths = find_paired_and_unpaired_files(fq_dir)
    stats = Counter()
    lengths = Counter()
    tile_nums = set()

    # Paired end reads
    ainr1, ainr2 = build_adapters()
    a_pre_r1 = dna_rev_comp(ainr2)
    a_pre_r2 = dna_rev_comp(ainr1)

    # Restriction site regex
    if rest_site == dna_rev_comp(rest_site):  # if palendromic
        pattern = rest_site
    else:
        pattern = '(%s|%s)' % (rest_site, dna_rev_comp(rest_site))
    rest_re = re.compile(pattern)

    def PE_restriction_sites(read_obj):
        if isinstance(read_obj, str):
            return [m.start() for m in rest_re.finditer(a_pre_r1 + read_obj + ainr1)]
        else:
            assert len(read_obj) == 2, read_obj
            read_obj = [a_pre_r1 + read_obj[0], a_pre_r2 + read_obj[1]]
            read1_pos = [m.start() for m in rest_re.finditer(read_obj[0])] 
            read2_pos = [len(read_obj[0]) + len(read_obj[1])  # total len
                            - m.start() - len(rest_site)  # minus that known to be to the right
                            for m in rest_re.finditer(read_obj[1])]
            return read1_pos + read2_pos

    def SE_restriction_sites(read):
        return [m.start() for m in rest_re.finditer(a_pre_r1 + read)]

    for fpath1, fpath2 in pe_fpaths:
        print os.path.basename(fpath1), os.path.basename(fpath2)
        it = itertools.izip(SeqIO.parse(gzip_friendly_open(fpath1), 'fastq'),
                            SeqIO.parse(gzip_friendly_open(fpath2), 'fastq'))
        if max_reads_per_file is not None:
            it = itertools.islice(it, None, max_reads_per_file)
        for r1, r2 in it:
            stats['total reads'] += 1
            assert r1.id == r2.id, '%s\n%s' % (r1.id, r2.id)
            tile_nums.add(int(str(r1.id).split(':')[4]))

            read_obj = collapse_if_overlapped_pair([str(r1.seq), str(r2.seq)])
            if isinstance(read_obj, str):
                lengths[len(read_obj)] += 1
            else:
                lengths['>%d' % sum(len(s) for s in read_obj)] += 1

            positions = PE_restriction_sites(read_obj)
            if positions:
                stats['site present'] += 1
            else:
                stats['site absent'] += 1

    for fpath in se_fpaths:
        print os.path.basename(fpath)
        it = SeqIO.parse(gzip_friendly_open(fpath), 'fastq')
        if max_reads_per_file is not None:
            it = itertools.islice(it, None, max_reads_per_file)
        for rec in it:
            stats['total reads'] += 1
            tile_nums.add(int(str(rec.id).split(':')[4]))
            positions = SE_restriction_sites({1: str(rec.seq)})
            if positions:
                stats['site present'] += 1
            else:
                stats['site absent'] += 1

    if tile_nums <= set(range(1101, 1115) + range(2101, 2115)):
        version = 'v1/v2'
    elif tile_nums <= set(range(1101, 1120) + range(2101, 2120)):
        version = 'v3'
    else:
        version = 'Unknown:', tile_nums

    print
    print project_name
    print 'MiSeq Version:', version
    print
    print 'Most common lengths'
    print 'Length: Count'
    for l, count in sorted(lengths.items(), key=lambda tup: tup[1], reverse=True)[:10]:
        print l, count
    print
    print 'Fraction with %s site: %d / %d = %.3f' % (
        rest_site,
        stats['site present'],
        stats['total reads'],
        float(stats['site present']) / stats['total reads']
    )


if __name__ == '__main__':
    usage_fmt = '%s <project_name> <restriction_enzyme_name> <restriction_site>' % sys.argv[0]
    if len(sys.argv) != len(usage_fmt.split()):
        sys.exit('Usage: ' + usage_fmt)
    project_name = sys.argv[1]
    rest_name = sys.argv[2]
    rest_site = sys.argv[3].upper()
    count_restriction_sites(project_name, rest_name, rest_site)
