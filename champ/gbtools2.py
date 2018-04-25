import h5py
import pysam
import cachetools
import lomp
import progressbar
from champ.kd import fit_one_group_kd
MINIMUM_REQUIRED_COUNTS = 5


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
            yield pileup_column.pos, frozenset(query_names)


def determine_kds_of_reads(contig_pileup_data, concentrations, delta_y, read_name_intensities):
    contig, pileup_data = contig_pileup_data
    position_kds = {}
    # since there are often long stretches where we have the reads, we cache the most recently-used ones and
    # try to look up the KD we calculated before.
    cache = cachetools.LRUCache(maxsize=20)
    for position, query_names in pileup_data:
        result = cache.get(query_names)
        if result is None:
            intensities = []
            for name in query_names:
                intensity_gradient = read_name_intensities.get(name)
                if intensity_gradient is None:
                    continue
                intensities.append(intensity_gradient)
            if len(intensities) < MINIMUM_REQUIRED_COUNTS:
                continue
            try:
                result = fit_one_group_kd(intensities, concentrations, delta_y=delta_y)
            except Exception as e:
                continue
            if result is None:
                continue
            cache[query_names] = result
        position_kds[position] = result
    return contig, position_kds


def calculate_genomic_kds(bamfile, read_name_intensities_hdf5_filename, concentrations, delta_y):
    print("loading read name intensities")
    read_name_intensities = load_read_name_intensities(read_name_intensities_hdf5_filename)
    with pysam.Samfile(bamfile) as samfile:
        contigs = list(reversed(sorted(samfile.references)))[:100]

    pileup_data = {}
    print("loading pileup data")
    with progressbar.ProgressBar(max_value=len(contigs)) as pbar:
        for contig in pbar(contigs):
            pileup_data[contig] = list(iterate_pileups(bamfile, contig))

    print("calculating genomic kds")
    contig_position_kds = {}
    with progressbar.ProgressBar(max_value=len(pileup_data)) as pbar:
        for contig, data in pbar(lomp.parallel_map(pileup_data.items(),
                                                   determine_kds_of_reads,
                                                   args=(concentrations, delta_y, read_name_intensities),
                                                   process_count=4)):
            contig_position_kds[contig] = data
    return contig_position_kds
