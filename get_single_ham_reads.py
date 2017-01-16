import local_config
import sys
import seqtools


def get_perfect_target_reads(target, reads_by_seq_fpath, out_fpath):
    single_ham_seqs = seqtools.get_mismatch_seqs(target, 1)
    single_ham_reads = {seq: [] for seq in single_ham_seqs}
    for line in open(reads_by_seq_fpath):
        words = line.strip().split()
        seq = words[0]
        read_names = words[1:]
        for sinham_seq in single_ham_seqs:
            if sinham_seq in seq:
                single_ham_reads[sinham_seq].extend(read_names)
    with open(out_fpath, 'w') as out:
        for sinham_seq in sorted(single_ham_seqs):
            out.write('{}\t{}\n'.format(sinham_seq, '\t'.join(single_ham_reads[sinham_seq])))



if __name__ == '__main__':
    usg_fmt = '{} <target_name> <reads_by_seq_fpath> <out_fpath>'.format(sys.argv[0])
    if len(sys.argv) != len(usg_fmt.split()):
        sys.exit(usg_fmt)
    
    target_name, reads_by_seq_fpath, out_fpath = sys.argv[1:]
    get_perfect_target_reads(local_config.targets[target_name],
                     reads_by_seq_fpath,
                     out_fpath)
