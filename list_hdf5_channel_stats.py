import h5py
import sys

if __name__ == '__main__':
    usg_fmt = '{} <hdf5_file>'.format(sys.argv[0])
    if len(sys.argv) != len(usg_fmt.split()):
        sys.exit('Usage: {}'.format(usg_fmt))

    h5_fpath = sys.argv[1]

    with h5py.File(h5_fpath) as f:
        channels = f.keys()
        channels.sort()

        print h5_fpath
        for channel in channels:
            print '{}: {} images'.format(channel, len(f[channel].keys()))
