import sys
import os
from glob import glob
import itertools
from collections import defaultdict, Counter
from Bio import SeqIO
import pysam
from uncertainties import ufloat
from pysam import Samfile
import numpy as np
import time
import pickle
import h5py
import lomp


float_vector_dt = h5py.special_dtype(vlen=np.float)
int_vector_dt = h5py.special_dtype(vlen=np.int32)
gene_affinity_dt = np.dtype([('kds', float_vector_dt),
                             ('kd_high_errors', float_vector_dt),
                             ('kd_low_errors', float_vector_dt),
                             ('counts', int_vector_dt),
                             ('breaks', int_vector_dt)])


class GeneAffinity(object):
    def __init__(self, name, affinity_data, gene_start, gene_end, cds_parts):
        self._name = name
        self._kds = []
        self._kd_error = []  # standard error of the fit?
        self._counts = []
        self._start = min(gene_start, gene_end)
        self._end = max(gene_end, gene_end)
        self._exon_kds = None
        self._exon_kd_errors = None
        self._exon_counts = None
        self._breaks = []
        self._is_valid = False
        last_good_position = None
        coverage_counter = 0
        for position in range(self._start, self._end):
            position_data = affinity_data.get(position)
            if position_data is None:
                self._kds.append(None)
                self._kd_error.append(None)
                self._counts.append(0)
                if last_good_position is not None and last_good_position == position - 1:
                    # we just transitioned to a gap in coverage from a region with coverage
                    self._breaks.append(coverage_counter)
            else:
                self._is_valid = True
                kd, minus_err, plus_err, count = position_data
                self._kds.append(kd)
                self._kd_error.append((plus_err, minus_err))
                self._counts.append(count)
                last_good_position = position
                coverage_counter += 1
        self._cds_parts = cds_parts

    @property
    def is_valid(self):
        return self._is_valid

    @property
    def name(self):
        return self._name

    @property
    def breaks(self):
        """
        The rightmost position in a region where we have coverage. Typically these are the 3' boundary of an exon
        but our exon library also has some coverage of introns.

        """
        return self._breaks

    @property
    def counts(self):
        return self._counts

    @property
    def compressed_counts(self):
        return [count for count in self._counts if count > 0]

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def all_kds(self):
        """
        KDs at each position in the gene, including introns. Positions without data will be set to None.

        """
        return self._kds

    @property
    def all_kd_errors(self):
        return self._kd_error

    @property
    def all_positions(self):
        """
        Coordinates of each base pair in the gene's chromosome.

        """
        return np.arange(self._start, self._end)

    @property
    def compressed_kds(self):
        """
        All KDs for this gene, in order, regardless of whether or not they are in an intron. Essentially just
        removing all the gaps in our coverage for easier interpretation.

        """
        return [kd for kd in self.all_kds if kd is not None]

    @property
    def compressed_95_percent_ci(self):
        return [err for err in self.all_kd_errors if err is not None]

    @property
    def exon_kds(self):
        """
        KDs at each exonic position in the gene. Positions without data will be set to None.

        """
        if self._exon_kds is None:
            self._calculate_exon_data()
        return self._exon_kds

    @property
    def exon_counts(self):
        if self._exon_counts is None:
            self._calculate_exon_data()
        return self._exon_counts

    @property
    def exon_kd_errors(self):
        if self._exon_kd_errors is None:
            self._calculate_exon_data()
        return self._exon_kd_errors

    def _calculate_exon_data(self):
        exonic_positions = set()
        for start, stop in self._cds_parts:
            for position in range(start, stop):
                exonic_positions.add(position)
        self._exon_kds = []
        self._exon_kd_errors = []
        self._exon_counts = []
        for kd, error, position, count in zip(self.all_kds, self.all_kd_errors, self.all_positions, self.counts):
            if position in exonic_positions:
                self._exon_kds.append(kd)
                self._exon_kd_errors.append(error)
                self._exon_counts.append(count)

    @property
    def highest_affinity(self):
        return np.min([k for k in self._kds if k is not None])


def get_qualifier_force_single(feature, qual_name, allow_zero=False):
    try:
        qual = feature.qualifiers[qual_name]
    except KeyError:
        if allow_zero:
            return ''
        else:
            raise
    if len(qual) > 1:
        print("WARNING:", feature, qual_name, qual)
    return qual[0]


class GenBankCDS(object):
    def __init__(self, rec, cds_feature):
        self.gene_id = get_qualifier_force_single(cds_feature, 'gene')
        self.chrm = rec.id
        self.strand = cds_feature.location.strand
        assert self.strand in [-1, 1]
        self.get_parts_and_boundaries(cds_feature)
        self.length = sum(abs(part[0] - part[1]) + 1 for part in self.parts)
        self.codon_start = get_qualifier_force_single(cds_feature, 'codon_start')
        try:
            self.protein_id = get_qualifier_force_single(cds_feature, 'protein_id')
            self.translation = get_qualifier_force_single(cds_feature, 'translation')
        except KeyError:
            assert 'exception' in cds_feature.qualifiers, cds
            self.protein_id = None
            self.translation = None

    def get_parts_and_boundaries(self, cds_feature):
        self.parts = []
        self.boundaries = set()
        for loc in cds_feature.location.parts:
            self.parts.append((loc.nofuzzy_start, loc.nofuzzy_end, loc.strand))
            self.boundaries.update((loc.nofuzzy_start, loc.nofuzzy_end))

    def __str__(self):
        return '%s: length %d, %s' % (self.gene_id, self.length, self.parts)


class GenBankGene(object):
    def __init__(self, rec, gene_feature):
        self.gene_id = get_qualifier_force_single(gene_feature, 'gene')
        syn_str = get_qualifier_force_single(gene_feature, 'gene_synonym', allow_zero=True)
        syn_str.replace(' ', '')
        self.gene_synonyms = set(syn_str.split(';'))
        note = get_qualifier_force_single(gene_feature, 'note', allow_zero=True)
        self.is_readthrough = bool('readthrough' in note.lower())

        self.chrm = rec.id
        self.gene_start = gene_feature.location.nofuzzy_start
        self.gene_end = gene_feature.location.nofuzzy_end

        self.cdss = []
        self.cds_parts = set()
        self.cds_boundaries = set()
        self.longest_cds = None
        self.longest_cds_length = None

    def contains_cds(self, cds):
        for part in cds.parts:
            if not (self.gene_start <= part[0] <= self.gene_end
                    and self.gene_start <= part[1] <= self.gene_end):
                return False
        return True

    def add_cds(self, cds):
        # First assert that gene id is the same and all parts are contained in the gene.
        assert cds.gene_id == self.gene_id or cds.gene_id in self.gene_synonyms, '%s\n%s' % (self, cds)
        assert self.contains_cds(cds), '%s\n%s' % (self, cds)

        self.cdss.append(cds)
        self.cds_parts.update(cds.parts)
        self.cds_boundaries.update(cds.boundaries)
        if self.longest_cds_length is None or cds.length > self.longest_cds_length:
            self.longest_cds = cds
            self.longest_cds_length = cds.length

    def overlaps_region(self, chrm, start, end):
        if self.chrm == chrm \
                and (self.gene_start <= start <= self.gene_end
                     or self.gene_start <= end <= self.gene_end
                     or start <= self.gene_start <= end
                     or start <= self.gene_end <= end):
            return True
        else:
            return False

    def overlaps_gene(self, other_gene):
        return self.overlaps_region(other_gene.chrm, other_gene.gene_start, other_gene.gene_end)

    def __str__(self):
        if self.longest_cds:
            longest_prot_id = self.longest_cds.protein_id
            all_prot_ids = ','.join([cds.protein_id for cds in self.cdss])
        else:
            longest_prot_id = None
            all_prot_ids = None
        return '%s\t%s:%d-%d\tnum_CDSs=%s\tlongest_CDS_protein=%s\tlongest_CDS_len=%s\tall_CDS_proteins=%s' \
               % (self.gene_id,
                  self.chrm,
                  self.gene_start,
                  self.gene_end,
                  len(self.cdss),
                  longest_prot_id,
                  self.longest_cds_length,
                  all_prot_ids)


def parse_gbff(fpath):
    """
    This method reads through refseq *_genomic.gbff files and extracts CDS isoform information.
    """
    assert fpath.endswith('.gbff'), 'Incorrect file type.'
    # Read in all the genes
    print('Reading genes from genomic file')
    genes_given_id = defaultdict(list)
    readthrough_genes = set()
    for rec in SeqIO.parse(open(fpath), 'gb'):
        if ('FIX_PATCH' in rec.annotations['keywords']
                or 'NOVEL_PATCH' in rec.annotations['keywords']
                or 'ALTERNATE_LOCUS' in rec.annotations['keywords']):
            continue
        for feature in rec.features:
            if feature.type == 'gene':
                gene = GenBankGene(rec, feature)
                if gene.is_readthrough:
                    # We reject any readthrough genes.
                    readthrough_genes.add(gene.gene_id)
                else:
                    genes_given_id[gene.gene_id].append(gene)
            elif feature.type == 'CDS':
                cds = GenBankCDS(rec, feature)
                if cds.protein_id is None:
                    # Some CDSs have exceptions indicating no protein is directly coded. Skip those
                    continue
                for gene in genes_given_id[cds.gene_id]:
                    if gene.contains_cds(cds):
                        gene.add_cds(cds)
                        break
                else:
                    if cds.gene_id not in readthrough_genes:
                        sys.exit('Error: Did not find gene for cds:\n\t%s' % cds)

    # In my observations, all large gene groups which use the same exons should be called isoforms
    # of the same gene.  They appear to simply be the invention of over-ambitious scientists,
    # potentially for more genes by their name or something like that.
    # Find groups which share cds 'exons'
    genes_given_part = defaultdict(list)
    for gene in itertools.chain(*genes_given_id.values()):
        for part in gene.cds_parts:
            genes_given_part[part].append(gene)
    overlapping_proteins = set()
    for genes in genes_given_part.values():
        if len(genes) > 1:
            overlapping_proteins.add(', '.join(sorted([gene.gene_id for gene in genes])))
    # Keep only the longest member of the family.
    for family_str in overlapping_proteins:
        family = family_str.split(', ')
        if len(family) <= 2:
            continue
        family_gene_iter = itertools.chain(*[genes_given_id[gene_id] for gene_id in family])
        max_len_gene = max(family_gene_iter, key=lambda gene: gene.longest_cds_length)
        for gene_id in family:
            if gene_id != max_len_gene.gene_id:
                del genes_given_id[gene_id]

    # Remove genes with no CDSs
    gene_ids_to_delete = set()
    for gene_id in genes_given_id.keys():
        genes = genes_given_id[gene_id]
        genes_to_keep = []
        for i in range(len(genes)):
            if genes[i].longest_cds is not None:
                genes_to_keep.append(i)
        if not genes_to_keep:
            gene_ids_to_delete.add(gene_id)
        else:
            genes_given_id[gene_id] = [genes[i] for i in genes_to_keep]
    for gene_id in gene_ids_to_delete:
        del genes_given_id[gene_id]

    # Search for duplicates
    dist_num_genes_for_gene_id = Counter()
    genes_with_dups = set()
    for gene_id, genes in genes_given_id.items():
        dist_num_genes_for_gene_id[len(genes)] += 1
        if len(genes) > 1:
            genes_with_dups.add(gene_id)
    print('Distribution of number of genes per gene_id:')
    for k, v in dist_num_genes_for_gene_id.items():
        print('\t%s: %s' % (k, v))
    print('Genes with duplicates:', ', '.join(sorted(genes_with_dups)))

    for n, (name, gene) in enumerate(genes_given_id.items()):
        if not gene:
            continue
        gene = gene[0]
        yield name, gene.chrm, gene.gene_start, gene.gene_end, gene.cds_parts


def convert_gbff_to_hdf5(hdf5_filename=None, gbff_filename=None):
    """ This finds the positions of all coding sequences for all genes in a Refseq-formatted GBFF file and saves the
    results to HDF5. It should be run once per flow cell. """

    # Use default paths if none are given
    if hdf5_filename is None:
        hdf5_filename = os.path.join(os.path.expanduser('~'), '.local', 'champ', 'gene-boundaries.h5')
    if gbff_filename is None:
        gbff_filename = os.path.join(os.path.expanduser("~"), '.local', 'champ', 'human-genome.gbff')


    string_dt = h5py.special_dtype(vlen=str)
    bounds_dt = np.dtype([('gene_id', np.uint32),
                          ('name', string_dt),
                          ('contig', string_dt),
                          ('gene_start', np.uint64),
                          ('gene_end', np.uint64)])
    cds_parts_dt = np.dtype([('gene_id', np.uint32),
                             ('cds_start', np.uint64),
                             ('cds_stop', np.uint64)])

    bounds = []
    all_cds_parts = []
    with h5py.File(hdf5_filename, 'w') as h5:
        for n, (name, contig, gene_start, gene_end, cds_parts) in enumerate(parse_gbff(gbff_filename)):
            bounds.append((n, name, contig, gene_start, gene_end))
            for start, stop, _ in cds_parts:
                all_cds_parts.append((n, start, stop))

        bounds_dataset = h5.create_dataset('/bounds', (len(bounds),), dtype=bounds_dt)
        bounds_dataset[...] = bounds

        cds_parts_dataset = h5.create_dataset('/cds-parts', (len(all_cds_parts),), dtype=cds_parts_dt)
        cds_parts_dataset[...] = all_cds_parts


def load_gene_positions(hdf5_filename=None):
    if hdf5_filename is None:
        hdf5_filename = os.path.join(os.path.expanduser('~'), '.local', 'champ', 'gene-boundaries.h5')
    with h5py.File(hdf5_filename, 'r') as h5:
        cds_parts = defaultdict(list)
        for gene_id, cds_start, cds_stop in h5['cds-parts'][:]:
            cds_parts[gene_id].append((cds_start, cds_stop))
        for gene_id, name, contig, gene_start, gene_stop in h5['bounds'][:]:
            yield gene_id, name, contig, gene_start, gene_stop, cds_parts[gene_id]


def parse_gene_affinities(contig, gene_start, gene_stop, position_kds):
    affinity_data = position_kds[contig]
    coverage_counter = 0
    last_good_position = None
    kds = []
    kd_high_errors = []
    kd_low_errors = []
    counts = []
    breaks = []
    for position in range(gene_start, gene_stop):
        position_data = affinity_data.get(position)
        if position_data is None:
            kds.append(None)
            kd_high_errors.append(None)
            kd_low_errors.append(None)
            counts.append(0)
            if last_good_position is not None and last_good_position == position - 1:
                # we just transitioned to a gap in coverage from a region with coverage
                breaks.append(coverage_counter)
        else:
            kd, minus_err, plus_err, count = position_data
            kds.append(kd)
            kd_high_errors.append(plus_err)
            kd_low_errors.append(minus_err)
            counts.append(count)
            last_good_position = position
            coverage_counter += 1
    return kds, kd_high_errors, kd_low_errors, counts, breaks



def build_gene_affinities(genes, position_kds):
    for gene_id, name, contig, gene_start, gene_stop, cds_parts in genes:
        print(gene_id)
        kds, kd_high_errors, kd_low_errors, counts, breaks = parse_gene_affinities(contig,
                                                                                   gene_start,
                                                                                   gene_stop,
                                                                                   position_kds)
        yield (gene_id,
               (kds, kd_high_errors, kd_low_errors, counts, breaks))
                  # (np.array(tuple(kds), dtype=float_vector_dt),
                  #  np.array(tuple(kd_high_errors), dtype=float_vector_dt),
                  #  np.array(tuple(kd_low_errors), dtype=float_vector_dt),
                  #  np.array(tuple(counts), dtype=int_vector_dt),
                  #  np.array(tuple(breaks), dtype=int_vector_dt))


def save_gene_affinities(gene_affinities, hdf5_filename=None):
    hdf5_filename = hdf5_filename if hdf5_filename is not None else os.path.join('results', 'gene-affinities.h5')
    with h5py.File(hdf5_filename, 'w') as h5:
        group = h5.create_group('gene-affinities')
        for gene_id, (kds, kd_high_errors, kd_low_errors, counts, breaks) in gene_affinities:
            dataset = group.create_dataset(str(gene_id), (), dtype=gene_affinity_dt)
            dataset['kds'] = kds
            dataset['kd_high_errors'] = kd_high_errors
            dataset['kd_low_errors'] = kd_low_errors
            dataset['counts'] = counts
            dataset['breaks'] = breaks


# def load_gene_affinities(gene_affinities_hdf5_path, gene_boundaries_hdf5_filename=None):
#     if gene_boundaries_hdf5_filename is None:
#         gene_boundaries_hdf5_filename = os.path.join(os.path.expanduser('~'),
#                                                      '.local',
#                                                      'champ',
#                                                      'gene-boundaries.h5')
#     with h5py.File(gene_affinities_hdf5_path, 'r') as gaffh5, h5py.File(gene_boundaries_hdf5_filename, 'r') as boundh5:




def load_contig_names(h5):
    contig_names = {}
    data = h5['contigs'][:]
    for contig_id, name in data:
        contig_names[contig_id] = name
    return contig_names


def load_position_kds(h5, contig_names):
    kds = defaultdict(dict)
    data = h5['genome-kds'][:]
    for contig_id, position, kd, minus_err, plus_err, count in data:
        contig_name = contig_names[contig_id]
        kds[contig_name][position] = kd, minus_err, plus_err, count
    return kds


def main_gaff(genome_kds_hdf5_path):
    """ We assume that convert_gbff_to_hdf5() has already been run using the default file paths. """
    print("\n\nLOADING GENE POSITIONS")
    genes = load_gene_positions()
    print("\n\nLOADING CONTIGS AND KDS")
    with h5py.File(genome_kds_hdf5_path, 'r') as h5:
        # TODO: The code that writes these is in genome-kd-march20. That needs to be moved here.
        contig_names = load_contig_names(h5)
        position_kds = load_position_kds(h5, contig_names)
    print("\n\nBUILDING GENE AFFINITIES")
    gene_affinities = build_gene_affinities(genes, position_kds)
    print("\n\nSAVING GENE AFFINITIES")
    save_gene_affinities(gene_affinities)