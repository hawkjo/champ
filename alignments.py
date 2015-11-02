from Bio import pairwise2
import itertools


class Error(object):
    bases = 'ACGTN'
    def __init__(self, ref_pos, **kwargs):
        if isinstance(self, Ins):
            assert int(ref_pos) == ref_pos - 0.5
        else:
            assert int(ref_pos) == ref_pos
            ref_pos = int(ref_pos)
        self.ref_pos = ref_pos


class SNP(Error):
    def __init__(self, ref_nt, new_nt, **kwargs):
        for n in [ref_nt, new_nt]:
            assert n in self.bases, n
        self.type = 'snp'
        self.ref_nt = ref_nt
        self.new_nt = new_nt
        super(SNP, self).__init__(**kwargs)

    def __str__(self):
        return 'SNP at %d: %s->%s' % (self.ref_pos, self.ref_nt, self.new_nt)


class Del(Error):
    def __init__(self, ref_nt, **kwargs):
        self.type = 'del'
        self.ref_nt = ref_nt
        super(Del, self).__init__(**kwargs)

    def __str__(self):
        return 'Del at %d: %s->\'-\'' % (self.ref_pos, self.ref_nt)


class Ins(Error):
    def __init__(self, new_nt, **kwargs):
        self.type = 'ins'
        self.new_nt = new_nt
        super(Ins, self).__init__(**kwargs)

    def __str__(self):
        return 'Ins at %.1f: \'-\'->%s' % (self.ref_pos, self.new_nt)


class Ambiguous(Error):
    def __init__(self, ref_nt, new_nt, **kwargs):
        self.type = 'amb'
        self.ref_nt = ref_nt
        self.new_nt = new_nt
        super(Ambiguous, self).__init__(**kwargs)

    def __str__(self):
        return 'Amb at %d: %s->%s' % (self.ref_pos, self.ref_nt, self.new_nt)


class Alignment:
    bases = 'ACGTN'

    submat = {}
    for c1, c2 in itertools.product(bases, bases):
        tup = (c1, c2)
        if 'N' in tup:
            submat[tup] = 0
        elif c1 == c2:
            submat[tup] = 2
        else:
            submat[tup] = -1

    def __init__(self, ref, read, read_name, max_indels=1):
        self.ref = ref
        self.read = read
        self.read_name = read_name
        self.max_indels = max_indels
        self._align()

    def _align(self):
        alignments = pairwise2.align.localds(self.ref, self.read, self.submat, -2, -1)
        # Force a single good alignment
        if len(alignments) == 1:
            self.alignment = alignments[0]
        else:
            self.alignment = None
            return

        self.full_ref_al, self.full_read_al, self.score, self.start, self.end = self.alignment
        self.ref_al = self.full_ref_al[self.start:self.end]
        self.read_al = self.full_read_al[self.start:self.end]
        if len(self.ref_al.replace('-', '')) == len(self.ref):
            self.end_to_end_ref_al = True

        self.errors = []
        self.Ns_in_read = []
        ref_pos = 0
        for refc, readc in zip(self.ref_al, self.read_al):
            if refc == 'N':
                if readc == '-':
                    self.errors.append(Del(ref_pos=ref_pos, ref_nt=refc))
                elif readc == 'N':
                    self.errors.append(Ambiguous(ref_pos=ref_pos, ref_nt=refc, new_nt=readc))
                self.Ns_in_read.append(readc)
            elif refc == readc:
                assert readc != '-', self.alignment
            elif refc == '-':
                ref_pos -= 1  # Undo advance of position
                self.errors.append(Ins(ref_pos=ref_pos + 0.5, new_nt=readc))
            elif readc == '-':
                self.errors.append(Del(ref_pos=ref_pos, ref_nt=refc))
            elif readc not in 'ACGT':
                self.errors.append(Ambiguous(ref_pos=ref_pos, ref_nt=refc, new_nt=readc))
            else:
                self.errors.append(SNP(ref_pos=ref_pos, ref_nt=refc, new_nt=readc))
            ref_pos += 1

    @property
    def is_perfect(self):
        return bool(len(self.errors) == 0 and len(self.ref_al.replace('-', '')) == len(self.ref))

    def print_errors(self):
        for error in self.errors:
            print error

    def __str__(self):
        return '\n'.join([self.ref, self.ref_al, self.read_al])
