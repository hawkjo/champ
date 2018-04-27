import os
from collections import defaultdict

import cachetools
import h5py
import lomp
from joblib import Parallel, delayed
import numpy as np
import progressbar
import pysam

from champ.kd import fit_one_group_kd

MINIMUM_REQUIRED_COUNTS = 5
QUALITY_THRESHOLD = 20
MAXIMUM_REALISTIC_DNA_LENGTH = 1000


def load_read_name_intensities(hdf5_filename):
    read_name_intensities = {}
    with h5py.File(hdf5_filename, 'r') as h5:
        read_names = h5['read_names'][:]
        intensity_matrix = h5['intensities'][:]
        for read_name, intensity_gradient in zip(read_names, intensity_matrix):
            read_name_intensities[read_name] = intensity_gradient
    return read_name_intensities


def iterate_pileups(bamfile, contig):
    """ Yields the sequence identifiers and the position they cover in the given contig. """
    with pysam.Samfile(bamfile) as sf:
        for pileup_column in sf.pileup(contig):
            if pileup_column.n < MINIMUM_REQUIRED_COUNTS:
                continue
            query_names = [pileup.alignment.query_name for pileup in pileup_column.pileups
                           if not pileup.alignment.is_qcfail and pileup.alignment.mapq > 20]
            yield contig, pileup_column.pos, frozenset(query_names)


def determine_kd_of_genomic_position(item, read_name_intensities, concentrations, delta_y):
    contig, position, query_names = item
    intensities = []
    for name in query_names:
        intensity_gradient = read_name_intensities.get(name)
        if intensity_gradient is None:
            continue
        intensities.append(intensity_gradient)
    if len(intensities) < MINIMUM_REQUIRED_COUNTS:
        return contig, position, None
    try:
        result = fit_one_group_kd(intensities, concentrations, delta_y=delta_y)
    except Exception:
        return contig, position, None
    return contig, position, result


# def determine_kds_of_reads(contig_pileup_data, concentrations, delta_y, read_name_intensities):
#     contig, pileup_data = contig_pileup_data
#     position_kds = {}
#     # since there are often long stretches where we have the reads, we cache the most recently-used ones and
#     # try to look up the KD we calculated before.
#     cache = cachetools.LRUCache(maxsize=20)
#     for position, query_names in pileup_data:
#         result = cache.get(query_names)
#         if result is None:
#             intensities = []
#             for name in query_names:
#                 intensity_gradient = read_name_intensities.get(name)
#                 if intensity_gradient is None:
#                     continue
#                 intensities.append(intensity_gradient)
#             if len(intensities) < MINIMUM_REQUIRED_COUNTS:
#                 continue
#             try:
#                 result = fit_one_group_kd(intensities, concentrations, delta_y=delta_y)
#             except Exception as e:
#                 continue
#             if result is None:
#                 continue
#             cache[query_names] = result
#         position_kds[position] = result
#     return contig, position_kds


def calculate_genomic_kds(bamfile, read_name_intensities_hdf5_filename, concentrations, delta_y):
    print("loading read name intensities")
    read_name_intensities = load_read_name_intensities(read_name_intensities_hdf5_filename)
    with pysam.Samfile(bamfile) as samfile:
        # contigs = list(reversed(sorted(samfile.references)))
        contigs = ['NC_000019.10']

    pileup_data = []
    print("loading pileup data")
    with progressbar.ProgressBar(max_value=len(contigs)) as pbar:
        for contig in pbar(contigs):
            for result in iterate_pileups(bamfile, contig):
                pileup_data.append(result)

    print("calculating genomic kds with lomp")
    contig_position_kds = {contig: {} for contig in contigs}
    with progressbar.ProgressBar(max_value=len(pileup_data)) as pbar:
        for contig, position, result in pbar(lomp.parallel_map(pileup_data,
                                                               determine_kd_of_genomic_position,
                                                               args=(read_name_intensities, concentrations, delta_y))):
            if result is not None:
                kd, kd_uncertainty, yint, fit_delta_y, count = result
                contig_position_kds[contig][position] = kd, kd_uncertainty, count
    return contig_position_kds


def load_gene_positions(hdf5_filename=None):
    if hdf5_filename is None:
        hdf5_filename = os.path.join(os.path.expanduser('~'), '.local', 'champ', 'gene-boundaries.h5')
    with h5py.File(hdf5_filename, 'r') as h5:
        cds_parts = defaultdict(list)
        for gene_id, cds_start, cds_stop in h5['cds-parts'][:]:
            cds_parts[gene_id].append((cds_start, cds_stop))
        for gene_id, name, sequence, contig, strand, gene_start, gene_stop in h5['bounds'][:]:
            yield gene_id, name, sequence, contig, strand, gene_start, gene_stop, cds_parts[gene_id]


def build_gene_affinities(genes, position_kds):
    for gene_id, name, sequence, contig, strand, gene_start, gene_stop, cds_parts in genes:
        if contig not in position_kds:
            continue
        kds, kd_uncertainties, counts, breaks = parse_gene_affinities(contig, gene_start, gene_stop, position_kds)
        yield (gene_id,
               np.array(kds, dtype=np.float),
               np.array(kd_uncertainties, dtype=np.float),
               np.array(counts, dtype=np.int32),
               np.array(breaks, dtype=np.int32))


def parse_gene_affinities(contig, gene_start, gene_stop, position_kds):
    affinity_data = position_kds[contig]
    coverage_counter = 0
    last_good_position = None
    kds = []
    kd_uncertainties = []
    counts = []
    breaks = []
    for position in range(gene_start, gene_stop):
        position_data = affinity_data.get(position)
        if position_data is None:
            kds.append(None)
            kd_uncertainties.append(None)
            counts.append(0)
            if last_good_position is not None and last_good_position == position - 1:
                # we just transitioned to a gap in coverage from a region with coverage
                breaks.append(coverage_counter)
        else:
            kd, kd_uncertainty, count = position_data
            kds.append(kd)
            kd_uncertainties.append(kd_uncertainty)
            counts.append(count)
            last_good_position = position
            coverage_counter += 1
    return kds, kd_uncertainties, counts, breaks


def load_gene_count(hdf5_filename=None):
    if hdf5_filename is None:
        hdf5_filename = os.path.join(os.path.expanduser('~'), '.local', 'champ', 'gene-boundaries.h5')
    with h5py.File(hdf5_filename, 'r') as h5:
        return len(h5['bounds'])


def save_gene_affinities(gene_affinities, gene_count, hdf5_filename=None):
    hdf5_filename = hdf5_filename if hdf5_filename is not None else os.path.join('results', 'gene-affinities.h5')
    with h5py.File(hdf5_filename, 'w') as h5:
        kd_dataset = h5.create_dataset('kds', (gene_count,),
                                       dtype=h5py.special_dtype(vlen=np.dtype('float32')))
        kd_uncertainties_dataset = h5.create_dataset('kd_uncertainties', (gene_count,),
                                                   dtype=h5py.special_dtype(vlen=np.dtype('float32')))
        counts_dataset = h5.create_dataset('counts', (gene_count,),
                                           dtype=h5py.special_dtype(vlen=np.dtype('int32')))
        breaks_dataset = h5.create_dataset('breaks', (gene_count,),
                                           dtype=h5py.special_dtype(vlen=np.dtype('int32')))

        for gene_id, kds, kd_high_errors, kd_low_errors, counts, breaks in gene_affinities:
            kd_dataset[gene_id] = kds
            kd_uncertainties_dataset[gene_id] = kd_high_errors
            counts_dataset[gene_id] = counts
            breaks_dataset[gene_id] = breaks


def main(bamfile, read_name_intensities_hdf5_filename, concentrations, delta_y,
         gene_boundaries_h5_path=None, gene_affinities_path=None):
    """ We assume that convert_gbff_to_hdf5() has already been run using the default file paths. """
    position_kds = calculate_genomic_kds(bamfile, read_name_intensities_hdf5_filename, concentrations, delta_y)
    genes = load_gene_positions(hdf5_filename=gene_boundaries_h5_path)
    gene_affinities = build_gene_affinities(genes, position_kds)
    gene_count = load_gene_count(hdf5_filename=gene_boundaries_h5_path)
    save_gene_affinities(gene_affinities, gene_count, hdf5_filename=gene_affinities_path)
