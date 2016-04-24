import h5py
import numpy as np
import os
import sys
import hdf5_tools
from subprocess import check_call
from misctools import cd


def sextractor_hdf5(h5_fpath):
    out_dir = os.path.splitext(h5_fpath)[0] + '_sex_files'
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    print 'Writing txt files...'
    all_out_fpaths = []
    with h5py.File(h5_fpath, 'r') as f:
        for channel, group in f.items():
            for dset_name, dset in group.items():
                sys.stdout.write('.')
                sys.stdout.flush()

                out_fname = hdf5_tools.bname_given_channel_and_dset_name(channel, dset_name) + '.txt'
                out_fpath = os.path.join(out_dir, out_fname)
                np.savetxt(out_fpath, dset, fmt='%d', delimiter='\t')
                all_out_fpaths.append(out_fpath)

    print '\nSextractor...'
    with cd(out_dir):
        for fpath in all_out_fpaths:
            bname = os.path.splitext(os.path.basename(fpath))[0]
            check_call(['spotproc.sh', bname])


if __name__ == '__main__':
    usage_fmt = '{} <h5_fpath>'.format(sys.argv[0])
    if len(sys.argv) != len(usage_fmt.split()):
        sys.exit(usage_fmt)
    sextractor_hdf5(sys.argv[1])
