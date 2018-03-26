import sys
import os
import lomp
from collections import defaultdict, Counter
from Bio import SeqIO
from pysam import Samfile
import numpy as np
import h5py
from scipy import stats
import progressbar
import itertools


MINIMUM_CLUSTER_COUNT = 6
QUALITY_THRESHOLD = 20
MAXIMUM_REALISTIC_DNA_LENGTH = 200


def load_genes_with_affinities(gene_boundaries_h5_filename=None, gene_affinities_filename=None):
    """ This is usually what you want to use when loading data. It will give you every gene with its
    KD data already loaded. """
    for gene in load_genes(gene_boundaries_h5_filename):
        kd, kd_low_errors, kd_high_errors, counts = load_gene_kd(gene.id, gene_affinities_filename)
        gene.set_measurements(kd, kd_low_errors, kd_high_errors, counts)
        yield gene


def load_genes(gene_boundaries_h5_filename=None):
    """ Normally you won't want to use this. It just loads genes and their positions (with exon data).
    This is something you'd only want to use if you wanted to examine the parsed Refseq data to see what's in
    that particular version. """
    if gene_boundaries_h5_filename is None:
        gene_boundaries_h5_filename = os.path.join(os.path.expanduser('~'), '.local', 'champ', 'gene-boundaries.h5')
    for gene_id, name, contig, gene_start, gene_stop, cds_parts in load_gene_positions(gene_boundaries_h5_filename):
        gaff = GeneAffinity(gene_id, name, contig)
        yield gaff.set_boundaries(gene_start, gene_stop, cds_parts)


def load_gene_kd(gene_id, hdf5_filename=None):
    """ Get the affinity data for a particular gene. """
    hdf5_filename = hdf5_filename if hdf5_filename is not None else os.path.join('results', 'gene-affinities.h5')
    with h5py.File(hdf5_filename, 'r') as h5:
        kd = h5['kds/{gene_id}'.format(gene_id=gene_id)].value
        kd_high_errors = h5['kd_high_errors/{gene_id}'.format(gene_id=gene_id)].value
        kd_low_errors = h5['kd_low_errors/{gene_id}'.format(gene_id=gene_id)].value
        counts = h5['counts/{gene_id}'.format(gene_id=gene_id)].value
        return kd, kd_low_errors, kd_high_errors, counts


class GeneAffinity(object):
    def __init__(self, gene_id, name, contig):
        """
        Represents our KD measurements across a gene. We might not have data at each location, especially if using
        an exome library.

        """
        self.id = gene_id
        self.name = name
        self.contig = contig
        self._kds = None
        self._kd_errors_low = None
        self._kd_errors_high = None
        self._counts = None
        self._exonic = None
        self._exon_boundaries = None

    def set_measurements(self, kds, kd_errors_low, kd_errors_high, counts):
        self._kds = kds
        self._kd_errors_low = kd_errors_low
        self._kd_errors_high = kd_errors_high
        self._counts = counts
        return self

    def set_boundaries(self, gene_start, gene_stop, cds_parts):
        gene_start, gene_stop = min(gene_start, gene_stop), max(gene_start, gene_stop)
        self._exonic = np.zeros(abs(gene_stop - gene_start), dtype=np.bool)
        self._exon_boundaries = []
        for cds_start, cds_stop in cds_parts:
            cds_start, cds_stop = min(cds_start, cds_stop), max(cds_start, cds_stop)
            start, stop = cds_start - gene_start, cds_stop - gene_start
            self._exonic[start:stop] = True
            self._exon_boundaries.append((start, stop))
        return self

    @property
    def exon_boundaries(self):
        for start, stop in self._exon_boundaries:
            yield start, stop

    def set_exons(self, exons):
        self._exonic = exons
        return self

    @property
    def kds(self):
        return self._kds

    @property
    def kd_errors_low(self):
        return self._kd_errors_low

    @property
    def kd_errors_high(self):
        return self._kd_errors_high

    @property
    def counts(self):
        return self._counts

    @property
    def exons(self):
        """ Collects data from exonic positions only, and repackages them into a Gene object. """
        positions = self._exonic
        kds = self._kds[positions]
        kd_errors_low = self._kd_errors_low[positions]
        kd_errors_high = self._kd_errors_high[positions]
        counts = self._counts[positions]
        exonic = np.ones(kds.shape, dtype=np.bool)
        gene = GeneAffinity(self.id, "%s Exons" % self.name, self.contig)
        gene = gene.set_measurements(kds, kd_errors_low, kd_errors_high, counts)
        gene = gene.set_exons(exonic)
        return gene

    @property
    def compressed(self):
        """ Collects data from positions where we made measurements (regardless of whether the position is exonic)
        and repackages them into a Gene object. """
        positions = ~np.isnan(self._kds)
        kds = self._kds[positions]
        kd_errors_low = self._kd_errors_low[positions]
        kd_errors_high = self._kd_errors_high[positions]
        counts = self._counts[positions]
        exonic = self._exonic[positions]
        gene = GeneAffinity(self.id, "%s Compressed" % self.name, self.contig)
        gene = gene.set_measurements(kds, kd_errors_low, kd_errors_high, counts)
        gene = gene.set_exons(exonic)
        return gene

    @property
    def highest_affinity(self):
        return np.nanmin(self._kds)


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
            assert 'exception' in cds_feature.qualifiers, cds_feature
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


def main_gaff(bamfile, read_name_kd_filename):
    """ We assume that convert_gbff_to_hdf5() has already been run using the default file paths. """
    read_name_kds = load_kds(read_name_kd_filename)
    position_kds = calculate_genomic_kds(bamfile, read_name_kds)
    genes = load_gene_positions()
    gene_affinities = build_gene_affinities(genes, position_kds)
    save_gene_affinities(gene_affinities)


def load_kds(filename):
    read_name_kds = {}
    with open(filename) as f:
        for line in f:
            read_name, kdstr, kderrstr = line.strip().split("\t")
            kd, kderr = float(kdstr), float(kderrstr)
            read_name_kds[read_name] = kd
    return read_name_kds


def calculate_genomic_kds(bamfile, read_name_kds):
    """

    """
    position_kds = {}
    try:
        with Samfile(bamfile) as samfile:
            # contigs = list(reversed(sorted(samfile.references)))
            contigs = ["NC_000019.10"]
            for n, contig in enumerate(contigs):
                contig_position_kds = find_kds_at_all_positions(samfile.fetch(contig), read_name_kds)
                position_kds[contig] = contig_position_kds
                print("%d/%d contigs complete." % (n, len(contigs)))
        return position_kds
    except IOError:
        raise ValueError("Could not open %s. Does it exist and is it valid?" % bamfile)


def load_gene_positions(hdf5_filename=None):
    if hdf5_filename is None:
        hdf5_filename = os.path.join(os.path.expanduser('~'), '.local', 'champ', 'gene-boundaries.h5')
    with h5py.File(hdf5_filename, 'r') as h5:
        cds_parts = defaultdict(list)
        for gene_id, cds_start, cds_stop in h5['cds-parts'][:]:
            cds_parts[gene_id].append((cds_start, cds_stop))
        for gene_id, name, contig, gene_start, gene_stop in h5['bounds'][:]:
            yield gene_id, name, contig, gene_start, gene_stop, cds_parts[gene_id]


def build_gene_affinities(genes, position_kds):
    for gene_id, name, contig, gene_start, gene_stop, cds_parts in genes:
        if contig not in position_kds:
            continue
        kds, kd_high_errors, kd_low_errors, counts, breaks = parse_gene_affinities(contig,
                                                                                   gene_start,
                                                                                   gene_stop,
                                                                                   position_kds)
        yield (gene_id,
               np.array(kds, dtype=np.float),
               np.array(kd_high_errors, dtype=np.float),
               np.array(kd_low_errors, dtype=np.float),
               np.array(counts, dtype=np.int32),
               np.array(breaks, dtype=np.int32))


def save_gene_affinities(gene_affinities, hdf5_filename=None):
    hdf5_filename = hdf5_filename if hdf5_filename is not None else os.path.join('results', 'gene-affinities.h5')
    with h5py.File(hdf5_filename, 'w') as h5:
        # But do that later
        kd_group = h5.create_group('kds')
        kd_high_errors_group = h5.create_group('kd_high_errors')
        kd_low_errors_group = h5.create_group('kd_low_errors')
        counts_group = h5.create_group('counts')
        breaks_group = h5.create_group('breaks')

        for gene_id, kds, kd_high_errors, kd_low_errors, counts, breaks in gene_affinities:
            gene_id_str = str(gene_id)
            kd_group.create_dataset(gene_id_str, data=kds)
            kd_high_errors_group.create_dataset(gene_id_str, data=kd_high_errors)
            kd_low_errors_group.create_dataset(gene_id_str, data=kd_low_errors)
            counts_group.create_dataset(gene_id_str, data=counts)
            breaks_group.create_dataset(gene_id_str, data=breaks)


def find_kds_at_all_positions(alignments, read_name_kds):
    """
    We want to know the KD at each base pair of the genomic DNA.

    """
    position_kds = defaultdict(list)
    sketchy = 0
    normal = 0
    short = 1
    for alignment in alignments:
        if alignment.is_reverse or alignment.is_qcfail:
            continue
        kd = read_name_kds.get(alignment.query_name)
        if kd is None:
            continue
        if alignment.is_paired and alignment.template_length > len(alignment.query_sequence):
            normal += 1
            start = alignment.reference_start
            end = start + alignment.template_length
        elif alignment.reference_length == len(alignment.query_sequence):
            short += 1
            start = alignment.reference_start
            end = start + alignment.reference_length
        else:
            # The read is unpaired and the alignment was sketchy, so we can't trust this read
            sketchy += 1
            continue
        if abs(end-start) > MAXIMUM_REALISTIC_DNA_LENGTH:
            continue
        # This is a good quality read and we can make valid claims about the affinity between start and end
        for position in range(start, end):
            position_kds[position].append((kd, start, end))
    print("sketchy %d, short %d, normal %d" % (sketchy, short, normal))
    final_results = {}
    pbar = progressbar.ProgressBar(max_value=len(position_kds))
    for position, median, ci_minus, ci_plus, count in pbar(lomp.parallel_map(position_kds.items(),
                                                                             _thread_find_best_offset_kd,
                                                                             process_count=12)):
        final_results[position] = median, ci_minus, ci_plus, count
    print("final results found: %d" % len(final_results))
    return final_results


def _thread_find_best_offset_kd(position_kd_data):
    position, kd_data = position_kd_data
    return find_best_offset_kd(position, kd_data)


def find_best_offset_kd(position, kd_data):
    # Don't worry about offsets since count is so low
    kds = [kd for kd, _, _ in kd_data]
    if len(kds) < 6:
        return position, None, None, None, len(kds)
    confidence95minus, confidence95plus = stats.mstats.median_cihs(kds)
    return position, np.median(kds), confidence95minus, confidence95plus, len(kds)

#
# def find_best_offset_kd(position, kd_data):
#     left_offset_kds = defaultdict(list)
#     right_offset_kds = defaultdict(list)
#     # We consider all reads that contain this particular base
#     # However, we also try looking at different windows of reads by varying how far to the left and right
#     # we allow the reads to start or end, respectively. This way, if our position is very close to a very
#     # high affinity site, we'll be able to throw out the reads that contain that other site and only consider
#     # ones that start after it (or end before it), and thus don't physically contain it. This gives us higher
#     # resolution
#     for kd, start, end in kd_data:
#         for offset in range(0, 100, 5):
#             # TODO: This can probably be more efficient if we precalculate the valid range
#             # TODO: But that's probably not a big deal, this shouldn't take long
#             if start >= (position - offset):
#                 left_offset_kds[offset].append(kd)
#             if end <= (position + offset):
#                 right_offset_kds[offset].append(kd)
#     if not left_offset_kds and not right_offset_kds:
#         # TODO: Potential issue: If we just have really long reads around this position, and no ends are close
#         # then we won't be able to get any data about it, even though we have a measurement that includes this
#         # position. This is reasonable, since our resolution is very poor in such cases, but we are throwing
#         # away data.
#         return position, None, None, None, 0
#     # We look at all the windows of reads and find the highest KD. This gives us the tightest affinity that
#     # can be plausibly linked to this particular position while excluding nearby high affinity sites
#     max_kd = 0.0
#     best_kds = None
#     for offset, kds in itertools.chain(left_offset_kds.items(), right_offset_kds.items()):
#         median = np.median(kds)
#         if median > max_kd:
#             max_kd = median
#             best_kds = kds, median  # cache the median so we don't have to recalculate it
#     if best_kds is None:
#         return position, None, None, None, 0
#     kds, median = best_kds
#     if len(kds) < 6:
#         return position, None, None, None, len(kds)
#     confidence95minus, confidence95plus = stats.mstats.median_cihs(kds)
#     return position, median, confidence95minus, confidence95plus, len(kds)
