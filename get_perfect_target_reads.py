import local_config
import sys


def get_perfect_target_reads(target, reads_by_seq_fpath, out_fpath):
    with open(out_fpath, 'w') as out:
        for line in open(reads_by_seq_fpath):
            words = line.strip().split()
            seq = words[0]
            read_names = words[1:]
            if target in seq:
                out.write('\n'.join(read_names) + '\n')


if __name__ == '__main__':
    usg_fmt = '{} <target_name> <reads_by_seq_fpath> <out_fpath>'.format(sys.argv[0])
    if len(sys.argv) != len(usg_fmt.split()):
        sys.exit(usg_fmt)
    
    target_name, reads_by_seq_fpath, out_fpath = sys.argv[1:]
    get_perfect_target_reads(local_config.targets[target_name],
                     reads_by_seq_fpath,
                     out_fpath)
