import sys
import pysam
import misc
from champ.kd import IAKdData
from scipy.optimize import curve_fit
from collections import defaultdict
import pyfasta
import pysam


class GenomicSequence(object):
    __slots__ = '_read_name', '_upstream', '_downstream', '_isize', 'fasta_start', 'reference_id', 'isize'

    def __init__(self, read_name):
        self._read_name = read_name
        self._upstream = None
        self._downstream = None
        self.isize = None
        self.reference_id = None
        self.fasta_start = None

    def add_sequence(self, sequence, reference_id, fasta_start, isize):
        assert len(sequence) < 4000
        if isize < 0:
            self._downstream = sequence
        else:
            self.isize = isize
            self._upstream = sequence
            self.reference_id = reference_id
            self.fasta_start = fasta_start

    def get_full_sequence(self, fasta):
        """ Gets the full sequence of the genomic DNA strand on the chip. If it's longer than 150 bp, we have to
        use sequence information from the assembled genome."""
        if self._upstream is not None and self._downstream is not None:
            inferred_start = self.fasta_start + len(self._upstream)
            inferred_end = self.fasta_start + self.isize - len(self._downstream)
            known_read_size = len(self._upstream) + len(self._downstream)
            if self.isize > known_read_size:
                # We weren't able to read the entire strand, so we infer the missing segment from the assembled genome
                return ("%s%s%s" % (self._upstream, fasta[self.reference_id][inferred_start:inferred_end], self._downstream)).upper(), self.isize
            else:
                # When DNA is shorter than the paired end read length, we'll have overlap, so we don't need to look
                # anything up.
                return "%s%s" % (self._upstream, self._downstream[known_read_size-self.isize:]), self.isize
        return None


def get_genomic_read_sequences(fasta_path, bamfile_path):
    fasta = pyfasta.Fasta(fasta_path, key_fn=lambda k: k.split()[0])
    sam = pysam.Samfile(bamfile_path)
    data = {}
    for read in sam:
        # ignore DNAs shorter than a PAM + target sequence and impossibly large sequences
        if abs(read.isize) < 24 or read.isize > 5000:
            continue
        if read.qname not in data:
            data[read.qname] = GenomicSequence(read.qname)
        data[read.qname].add_sequence(read.seq, read.reference_name, read.reference_start, read.isize)

    for read_name, gs in data.items():
        result = gs.get_full_sequence(fasta)
        if result is not None:
            seq, size = result
            yield read_name, seq, size


class ScoredRead(object):
    """A container class usable within DoublyLinkedScoreLists"""
    def __init__(self,
                 read_name,
                 start,
                 end,
                 concs,
                 scores):
        self.read_name = read_name
        self.start = start
        self.end = end
        self.concs = concs
        self.scores = scores
        self.prev = None
        self.next = None


class DoublyLinkedScoreList(object):
    """A doubly-linked list with ScoredReads as nodes. O(1) append/remove sorted list given sorted input."""
    head = None
    tail = None
    _min_end = None
    _min_end_nodes = set()
    _len = 0

    def append(self, *args):
        new_node = ScoredRead(*args)
        if self.head is None:
            self.head = self.tail = new_node
        else:
            new_node.prev = self.tail
            self.tail.next = new_node
            self.tail = new_node

        if self._min_end is None:
            self._min_end = new_node.end
            self._min_end_nodes = set([new_node])
        elif new_node.end == self._min_end:
            self._min_end_nodes.add(new_node)
        elif new_node.end < self._min_end:
            self._min_end = new_node.end
            self._min_end_nodes = set([new_node])
        self._len += 1

    def remove(self, some_node):
        if some_node is self.head and some_node is self.tail:
            self.head = None
            self.tail = None
        elif some_node is self.head:
            self.head = some_node.next
            self.head.prev = None
        elif some_node is self.tail:
            self.tail = some_node.prev
            self.tail.next = None
        else:
            some_node.prev.next = some_node.next
            some_node.next.prev = some_node.prev
        self._len -= 1

        if some_node.end == self._min_end:
            self._min_end_nodes -= set([some_node])
            if not self._min_end_nodes:
                self._update_min_end()

    def remove_current_min_end_reads(self):
        mn_end_list = list(self._min_end_nodes)[:]
        for rs in mn_end_list:
            self.remove(rs)

    def _update_min_end(self):
        if self.head is None:
            self._min_end = None
        else:
            self._min_end = self.head.end
            self._min_end_nodes = set([self.head])
            for nd in self:
                if nd.end < self._min_end:
                    self._min_end = nd.end
                    self._min_end_nodes = set([nd])
                elif nd.end == self._min_end:
                    self._min_end_nodes.add(nd)

    def __iter__(self):
        current_node = self.head
        while current_node is not None:
            yield current_node
            current_node = current_node.next

    def reverse_iter(self):
        current_node = self.tail
        while current_node is not None:
            yield current_node
            current_node = current_node.prev

    def __len__(self):
        return self._len

    @property
    def min_end(self):
        return self._min_end


class KdFitGenome(object):
    """Class to fit Kds at every position of of interest in a genome."""
    # Class runs through reads in a given bam file, piling reads up for Kd fitting.
    #
    # At every position where at least one read starts or stops and enough reads are piled up, the
    # Kd value using all reads is calculated, as well as Kds with desired limiting of overhangs in
    # each direction.
    def __init__(self,
                 int_scores,
                 h5_fpaths,
                 signal_channel,
                 IA_Kd_fpath,
                 directional_Kd_offsets=[],
                 min_clust=5,
                 mapq_cutoff=20):
        self.int_scores = int_scores
        self.h5_fpaths = h5_fpaths
        self.signal_channel = signal_channel

        self.all_read_names = set()
        for h5_fpath in self.h5_fpaths:
            self.all_read_names.update(
                self.int_scores.score_given_read_name_in_channel[h5_fpath][signal_channel].keys()
            )
        self.concentrations = map(misc.parse_concentration, self.h5_fpaths)

        self.IAKdData = IAKdData(IA_Kd_fpath)
        assert self.concentrations == self.IAKdData.concentrations, (self.concentrations,
                                                                   self.IAKdData.concentrations)
        self.Imin = self.IAKdData.Imin
        self.Imax = self.IAKdData.Imax
        assert len(self.Imin) == len(self.Imax), (self.Imin, self.Imax)
        self.Irange = [float(imx - imn) for imn, imx in zip(self.Imin, self.Imax)]

        self.directional_Kd_offsets = directional_Kd_offsets
        self.num_outputs_per_pos = 1 + 4 * len(self.directional_Kd_offsets)

        self.min_clust = min_clust
        self.mapq_cutoff = mapq_cutoff

    def add_read_scores_to_list(self, read_name, start, end):
        """Build ReadScore and add to list"""
        concs, read_scores = [], []
        for h5_fpath, conc, imn, irng in zip(self.h5_fpaths,
                                            self.concentrations,
                                            self.Imin,
                                            self.Irange):
            score_dict = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
            if read_name in score_dict:
                concs.append(conc)
                read_scores.append((score_dict[read_name] - imn)/irng)
        self.read_scores_list.append(read_name,
                                     start,
                                     end,
                                     concs,
                                     read_scores)

    def remove_passed_read_scores(self, pos):
        """Remove reads no longer relevant"""
        for rs in self.read_scores_list:
            if rs.end <= pos:
                self.read_scores_list.remove(rs)

    def Iobs(self, x, Kd):
        return 1.0/(1 + (float(Kd)/x))

    def fit_one_Kd(self, concs, scores):
        popt, pcov = curve_fit(self.Iobs, concs, scores, maxfev=100000)
        return popt[0]

    def fit_Kds_at_pos(self, pos, out_fh):
        """Fit Kds using all reads and all requested directional subsets overlapping pos"""
        # Housekeeping
        self.remove_passed_read_scores(pos)
        if len(self.read_scores_list) < self.min_clust:
            if self.last_write_contained_Kds:
                out_fh.write('{:d}\t'.format(pos) + '\t'.join('-' * self.num_outputs_per_pos) + '\n')
            self.last_write_contained_Kds = False
            return

        # Fit Kd with all
        all_concs = [conc for rs in self.read_scores_list for conc in rs.concs]
        all_scores = [score for rs in self.read_scores_list for score in rs.scores]
        outputs = [self.fit_one_Kd(all_concs, all_scores), len(all_scores)]

        # Fit "directional" Kds, filtering reads hanging over too far in one direction
        for offset in self.directional_Kd_offsets:
            left_read_count, right_read_count = 0, 0
            left_concs, left_scores, right_concs, right_scores = [], [], [], []
            for rs in self.read_scores_list:
                if rs.end <= pos + offset:
                    left_read_count += 1
                    left_concs.extend(rs.concs)
                    left_scores.extend(rs.scores)
                if rs.start >= pos - offset:
                    right_read_count += 1
                    right_concs.extend(rs.concs)
                    right_scores.extend(rs.scores)

            if left_read_count >= self.min_clust:
                outputs.append(self.fit_one_Kd(left_concs, left_scores))
                outputs.append(left_read_count)
            else:
                outputs.extend(['-', '-'])
            if right_read_count >= self.min_clust:
                outputs.append(self.fit_one_Kd(right_concs, right_scores))
                outputs.append(right_read_count)
            else:
                outputs.extend(['-', '-'])

        # Write results
        out_fh.write('{:d}\t'.format(pos) + '\t'.join('{}'.format(val) for val in outputs) + '\n')
        self.last_write_contained_Kds = True

    def finish_contig_Kds(self, start_pos, out_fh):
        """After reading in last read in contig, fit Kd at all remaining locations"""
        # Find remaining positions
        remaining_pos = set([start_pos])
        for rs in self.read_scores_list:
            for pos in [rs.start, rs.end]:
                if pos >= start_pos:
                    remaining_pos.add(pos)
        remaining_pos = list(sorted(remaining_pos))

        # Fit Kds at all remaining positions of interest
        for pos in remaining_pos:
            self.fit_Kds_at_pos(pos, out_fh)
        out_fh.write('{:d}\t'.format(remaining_pos[-1]) + '\t'.join('-' * self.num_outputs_per_pos) + '\n')

        # Clean out any remaining reads
        for rs in self.read_scores_list:
            self.read_scores_list.remove(rs)

    def fit_Kds_in_bam_and_write_results(self, bam_fpath, out_fpath):
        """Fit Kds at every status change in overlapping reads"""
        self.ColumnTitles = ['Pos', 'Kd_All', 'Cov']
        for offset in self.directional_Kd_offsets:
            self.ColumnTitles.append('Kd_<=+{:d}bp'.format(offset))
            self.ColumnTitles.append('Cov')
            self.ColumnTitles.append('Kd_>=-{:d}bp'.format(offset))
            self.ColumnTitles.append('Cov')

        def read_qc_and_ends(read):
            if (read.is_qcfail
                or read.mapq < self.mapq_cutoff
                or read.qname not in self.all_read_names):
                return None
            if read.is_paired and read.isize > len(read.seq):
                start = read.pos
                end = start + read.isize
                read_len = read.isize
            elif read.alen == len(read.seq):
                start = read.pos
                end = start + read.alen
                read_len = read.alen
            else:
                return None
            return start, end, read_len

        with open(out_fpath, 'w') as out:
            # Headers
            out.write('# ' + '\t'.join(self.ColumnTitles) + '\n')

            # Initialize
            sf = pysam.Samfile(bam_fpath)
            qc_res = None
            while qc_res is None:
                read = next(sf)
                qc_res = read_qc_and_ends(read)
            start, end, read_len = qc_res

            self.read_scores_list = DoublyLinkedScoreList()
            self.add_read_scores_to_list(read.qname, start, end)
            prev_start = start
            prev_contig = read.rname
            self.read_lens = [read_len]
            out.write('>{}\n'.format(prev_contig))
            self.last_write_contained_Kds = False

            # Proceed
            for i, read in enumerate(sf):
                if i % 10000 == 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
                qc_res = read_qc_and_ends(read)
                if qc_res is None:
                    continue
                start, end, read_len = qc_res
                self.read_lens.append(read_len)
                self.curr_contig = read.rname
                if self.curr_contig != prev_contig:
                    self.finish_contig_Kds(prev_start, out)
                    prev_contig = self.curr_contig
                    out.write('>{}\n'.format(prev_contig))
                    self.last_write_contained_Kds = False
                    sys.stdout.write('*')
                    sys.stdout.flush()
                elif start != prev_start:
                    self.fit_Kds_at_pos(prev_start, out)
                    while (self.read_scores_list.min_end is not None
                           and self.read_scores_list.min_end < start):
                        self.fit_Kds_at_pos(self.read_scores_list.min_end, out)
                        self.read_scores_list.remove_current_min_end_reads()
                    prev_start = start
                self.add_read_scores_to_list(read.qname, start, end)
            self.finish_contig_Kds(start, out)


class KdGenomeData(object):
    def __init__(self, Genome_Kd_fpath, IA_Kd_fpath):
        self.fpath = Genome_Kd_fpath
        self.IAKdData = IAKdData(IA_Kd_fpath)

    @property
    def all_full_Kds(self):
        for line in open(self.fpath):
            if not line.startswith('#') and not line.startswith('>'):
                words = line.strip().split()
                if not len(words) > 1:
                    continue
                assert len(words) > 1, words
                if words[1] != '-':
                    yield float(words[1])

    @property
    def all_full_ABAs(self):
        for Kd in self.all_full_Kds:
            yield self.IAKdData.ABA_given_Kd(Kd)

    def load_Kds(self):
        self.locs = defaultdict(list)
        self.Kds = defaultdict(list)
        self.coverage = defaultdict(list)
        self.max_Kds = defaultdict(list)
        self.max_Kd_coverage = defaultdict(list)
        for line in open(self.fpath):
            if line.startswith('#'):
                continue
            elif line.startswith('>'):
                curr_chrm = line.strip().split()[0][1:] # first word minus '>'
            else:
                words = line.strip().split()
                loc = int(words[0])
                all_Kds = [float(words[i]) if words[i] != '-' else None for i in range(1, len(words), 2)]
                all_covs = [int(words[i]) if words[i] != '-' else 0 for i in range(2, len(words), 2)]
                max_Kd_idx = max(range(len(all_Kds)), key=lambda i: all_Kds[i])
                max_Kd = all_Kds[max_Kd_idx]
                max_Kd_cov = all_covs[max_Kd_idx]

                self.locs[curr_chrm].append(loc)
                self.Kds[curr_chrm].append(all_Kds[0])
                self.coverage[curr_chrm].append(all_covs[0])
                self.max_Kds[curr_chrm].append(max_Kd)
                self.max_Kd_coverage[curr_chrm].append(max_Kd_cov)
        self.locs = dict(self.locs)
        self.Kds = dict(self.Kds)
        self.coverage = dict(self.coverage)
        self.max_Kds = dict(self.max_Kds)
        self.max_Kd_coverage = dict(self.max_Kd_coverage)

    def Kds_in_range(self, chrm, start, end, max_Kds=False):
        assert start < end
        if chrm not in self.locs:
            return [], [], []

        start_idx, end_idx = 0, 0
        for i, loc in enumerate(self.locs[chrm]):
            if loc <= start:
                start_idx = i
            elif loc > end:
                end_idx = i
                break
        if end_idx == 0:
            end_idx = -1
        if max_Kds:
            return (self.locs[chrm][start_idx:end_idx],
                    self.max_Kds[chrm][start_idx:end_idx],
                    self.max_Kd_coverage[chrm][start_idx:end_idx])
        else:
            return (self.locs[chrm][start_idx:end_idx],
                    self.Kds[chrm][start_idx:end_idx],
                    self.coverage[chrm][start_idx:end_idx])

    def ABAs_in_range(self, chrm, start, end, min_ABAs=False):
        locs, Kds, coverage = self.Kds_in_range(chrm, start, end, max_Kds=min_ABAs)
        ABAs = map(self.IAKdData.ABA_given_Kd, Kds)
        return locs, ABAs, coverage
