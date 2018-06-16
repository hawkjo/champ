import os
from collections import defaultdict

import h5py
import lomp
import numpy as np
import progressbar
import pysam
from more_itertools import windowed

from champ.constants import MINIMUM_REQUIRED_COUNTS
from champ.kd import determine_kd_of_genomic_position

QUALITY_THRESHOLD = 20
MAXIMUM_REALISTIC_DNA_LENGTH = 1000


class GeneAffinity(object):
    def __init__(self, gene_id, name, contig, sequence):
        """
        Represents our KD measurements across a gene. We might not have data at each location, especially if using
        an exome library.
        """
        self.id = gene_id
        self.name = name
        self.contig = contig
        self.sequence = sequence
        self._kds = None
        self._counts = None
        self._exonic = None
        self._exon_boundaries = None
        self._compressed_boundaries = None
        self._strand = None

    def set_measurements(self, kds, counts):
        self._kds = kds
        self._counts = counts
        return self

    def set_boundaries(self, strand, gene_start, gene_stop, cds_parts):
        assert gene_start < gene_stop
        self._strand = strand
        self._exonic = np.zeros(gene_stop - gene_start, dtype=np.bool)
        self._exon_boundaries = []
        min_start, max_stop = None, None
        for cds_start, cds_stop in cds_parts:
            cds_start, cds_stop = min(cds_start, cds_stop), max(cds_start, cds_stop)
            assert cds_start < cds_stop
            start, stop = cds_start - gene_start, cds_stop - gene_start
            assert start < stop
            self._exonic[start:stop] = True
            self._exon_boundaries.append((start, stop))
            min_start = min(start, min_start) if min_start is not None else start
            max_stop = max(stop, max_stop) if max_stop is not None else stop
        return self
    
    @property
    def exon_boundaries(self):
        for start, stop in self._exon_boundaries:
            yield start, stop

    @property
    def compressed_boundaries(self):
        return self._compressed_boundaries

    def set_exons(self, exons):
        self._exonic = exons
        return self

    @property
    def kds(self):
        return self._kds

    @property
    def counts(self):
        return self._counts

    @property
    def exons(self):
        """ Collects data from exonic positions only, and repackages them into a Gene object. """
        positions = self._exonic
        kds = self._kds[positions]
        counts = self._counts[positions]
        exonic = np.ones(kds.shape, dtype=np.bool)
        sequence = None if self.sequence is None else ''.join([self.sequence[n] for n, i in enumerate(positions) if i])
        gene = GeneAffinity(self.id, "%s Exons" % self.name, self.contig, sequence)
        gene = gene.set_measurements(kds, counts)
        gene = gene.set_exons(exonic)
        return gene

    @property
    def compressed(self):
        """ Collects data from positions where we made measurements (regardless of whether the position is exonic)
        and repackages them into a Gene object. """
        positions = self._counts >= MINIMUM_REQUIRED_COUNTS
        kds = self._kds[positions]
        counts = self._counts[positions]
        exonic = self._exonic[positions]
        compressed_boundaries = []
        position = 0
        for first, second in windowed(self._counts, 2, fillvalue=0):
            if second >= MINIMUM_REQUIRED_COUNTS > first:
                compressed_boundaries.append(position)
            if first >= MINIMUM_REQUIRED_COUNTS:
                position += 1
        sequence = None if self.sequence is None else ''.join([self.sequence[n] for n, i in enumerate(positions) if i])
        gene = GeneAffinity(self.id, "%s" % self.name, self.contig, sequence)
        gene = gene.set_measurements(kds, counts)
        gene = gene.set_exons(exonic)
        gene._compressed_boundaries = compressed_boundaries
        return gene

    @property
    def highest_affinity(self):
        return np.nanmin(self._kds)

    @property
    def highest_affinity_location(self):
        return np.nanargmin(self._kds)

    @property
    def coverage_count(self):
        """ Number of valid measurements """
        return self._kds[~np.isnan(self._kds)].shape[0]

    @property
    def coverage_ratio(self):
        return float(self.coverage_count) / self._kds.shape[0]


def load_genes_with_affinities(gene_boundaries_h5_filename=None, gene_affinities_filename=None):
    """ This is usually what you want to use when loading data. It will give you every gene with its
    KD data already loaded. """
    gene_affinities_filename = gene_affinities_filename if gene_affinities_filename is not None else os.path.join('results', 'cluster-data.h5')
    with h5py.File(gene_affinities_filename, 'r') as h5:
        for gene in load_genes(gene_boundaries_h5_filename):
            kd, counts = load_gene_kd(h5, gene.id)
            if len(kd) == 0:
                continue
            gene.set_measurements(kd, counts)
            yield gene


def load_gene_kd(h5, gene_id):
    """ Get the affinity data for a particular gene. """
    kd = h5['genome-kds'][gene_id]
    counts = h5['genome-counts'][gene_id]
    return kd, counts


def load_genes(gene_boundaries_h5_filename=None):
    """ Normally you won't want to use this. It just loads genes and their positions (with exon data).
    This is something you'd only want to use if you wanted to examine the parsed Refseq data to see what's in
    that particular version. """
    if gene_boundaries_h5_filename is None:
        gene_boundaries_h5_filename = os.path.join(os.path.expanduser('~'), '.local', 'champ', 'gene-boundaries.h5')
    for gene_id, name, sequence, contig, strand, gene_start, gene_stop, cds_parts in load_gene_positions(gene_boundaries_h5_filename):
        gaff = GeneAffinity(gene_id, name, contig, sequence)
        yield gaff.set_boundaries(strand, gene_start, gene_stop, cds_parts)


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


def calculate_genomic_kds(bamfile, read_name_intensities_hdf5_filename, concentrations, delta_y, process_count=16):
    print("calculate_genomic_kds")
    read_name_intensities = load_read_name_intensities(read_name_intensities_hdf5_filename)
    with pysam.Samfile(bamfile) as samfile:
        contigs = list(reversed(sorted(samfile.references)))
        # contigs = ['NC_000019.10']

    contig_position_kds = {contig: {} for contig in contigs}
    for contig in contigs:
        if not contig.startswith('NC'):
            # only examine DNA aligned to full chromosomes
            continue
        pileup_data = []
        for result in iterate_pileups(bamfile, contig):
            pileup_data.append(result)
        if len(pileup_data) > 500000:
            print(contig)
            pbar = progressbar.ProgressBar(max_value=len(pileup_data))
        else:
            print("small contig: %s" % contig)
            pbar = lambda x: x
        for position, result in pbar(lomp.parallel_map(pileup_data,
                                                       determine_kd_of_genomic_position,
                                                       args=(read_name_intensities, concentrations, delta_y),
                                                       process_count=process_count)):
            if result is not None:
                kd, kd_uncertainty, yint, fit_delta_y, count = result
                contig_position_kds[contig][position] = kd, count
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
        kds, counts, breaks = parse_gene_affinities(contig, gene_start, gene_stop, position_kds)
        yield (gene_id,
               np.array(kds, dtype=np.float),
               np.array(counts, dtype=np.int32),
               np.array(breaks, dtype=np.int32))


def parse_gene_affinities(contig, gene_start, gene_stop, position_kds):
    affinity_data = position_kds[contig]
    coverage_counter = 0
    last_good_position = None
    kds = []
    # kd_uncertainties = []
    counts = []
    breaks = []
    for position in range(gene_start, gene_stop):
        position_data = affinity_data.get(position)
        if position_data is None:
            kds.append(None)
            # kd_uncertainties.append(None)
            counts.append(0)
            if last_good_position is not None and last_good_position == position - 1:
                # we just transitioned to a gap in coverage from a region with coverage
                breaks.append(coverage_counter)
        else:
            kd, count = position_data
            kds.append(kd)
            counts.append(count)
            last_good_position = position
            coverage_counter += 1
    return kds, counts, breaks


def load_gene_count(hdf5_filename=None):
    if hdf5_filename is None:
        hdf5_filename = os.path.join(os.path.expanduser('~'), '.local', 'champ', 'gene-boundaries.h5')
    with h5py.File(hdf5_filename, 'r') as h5:
        return len(h5['bounds'])


def save_gene_affinities(gene_affinities, gene_count, hdf5_filename=None):
    hdf5_filename = hdf5_filename if hdf5_filename is not None else os.path.join('results', 'cluster-data.h5')
    with h5py.File(hdf5_filename, 'a') as h5:
        kd_dataset = h5.create_dataset('genome-kds', (gene_count,),
                                       dtype=h5py.special_dtype(vlen=np.dtype('float32')))
        counts_dataset = h5.create_dataset('genome-counts', (gene_count,),
                                           dtype=h5py.special_dtype(vlen=np.dtype('int32')))
        breaks_dataset = h5.create_dataset('genome-breaks', (gene_count,),
                                           dtype=h5py.special_dtype(vlen=np.dtype('int32')))

        for gene_id, kds, counts, breaks in gene_affinities:
            kd_dataset[gene_id] = kds
            counts_dataset[gene_id] = counts
            breaks_dataset[gene_id] = breaks


def genome_main(bamfile, read_name_intensities_hdf5_filename, concentrations, median_saturated_intensity, process_count=16):
    position_kds = calculate_genomic_kds(bamfile, read_name_intensities_hdf5_filename, concentrations, median_saturated_intensity, process_count=process_count)
    genes = load_gene_positions()
    gene_affinities = build_gene_affinities(genes, position_kds)
    gene_count = load_gene_count()
    save_gene_affinities(gene_affinities, gene_count, hdf5_filename=read_name_intensities_hdf5_filename)


def get_quality_alignment_start_and_end(alignment):
    if alignment.is_qcfail or alignment.mapq < 20:
        return None
    if alignment.is_proper_pair:
        # This read name appears twice in our data - once in the forward and once in the reverse direction
        if alignment.is_reverse:
            # For paired end reads, the forward and reverse reads are symmetric, so to avoid double counting the
            # bases between the two reads, we only count the bases once (with the forward read)
            return None
        if alignment.template_length <= alignment.reference_length or alignment.reference_length == 0:
            # the combined length of the first and second read is no greater than just one read - either the reads
            # perfectly overlapped or something weird is happening. We discard this read to be safe.
            return None
        start = alignment.reference_start
        end = start + alignment.template_length
    # elif alignment.reference_length == alignment.query_length and alignment.reference_length > 0:
    elif alignment.reference_length > 0:
        # This is an unpaired read, which is rare but not unheard of. We require that the read align perfectly to the
        # reference sequence, so we aren't dealing with indels. This might be a bit conservative - we'll come back to
        # this decision if the read counts are terrible
        start = alignment.reference_start
        end = start + alignment.reference_length
        assert start < end
    elif alignment.reference_length == 0:
        return None
    elif alignment.reference_length < 0:
        return None
    else:
        # The read is unpaired and the alignment was sketchy, so we can't trust this read
        return None
    if abs(end - start) > MAXIMUM_REALISTIC_DNA_LENGTH:
        return None
    return start, end
