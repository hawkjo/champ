import tifffile
import h5py
import hdf5_tools
import sys
import os
import glob
import re
import numpy as np


def tif_dir_to_hdf5(Major_axis_idx, h5_bname, tif_fpaths):
    Major_axis_group = Major_axis_idx + 1
    minor_axis_group = (1 - Major_axis_idx) + 1
    for fpath in tif_fpaths:
        sys.stdout.write('.')
        sys.stdout.flush()

        fname = os.path.split(fpath)[1]
        m = re.search('Pos_(\d+)_(\d+)\D', fname)
        Major_axis_pos = int(m.group(Major_axis_group))
        minor_axis_pos = int(m.group(minor_axis_group))

        with tifffile.TiffFile(fpath) as tif:
            summary = tif.micromanager_metadata['summary']

            # Find channel names and assert unique 
            channel_names = [name.replace(' ', '_').replace('(', '').replace(')', '')
                             for name in summary['ChNames']]
            assert summary['Channels'] == len(channel_names) == len(set(channel_names)), channel_names

            # channel_idxs map tif pages to channels
            channel_idxs = tif.micromanager_metadata['index_map']['channel']
            frame_idxs = tif.micromanager_metadata['frame']

            # Setup summed_images
            h, w = summary['Height'], summary['Width']
            def hw_zeros():
                return np.zeros((h, w), dtype=np.int)
            summed_images = {frame_idx: {channel_idx: hw_zeros() for channel_idx in channel_idxs}
                             for frame_idx in frame_idxs}

            # Add images
            array_lens = map(len, (channel_idxs, frame_idxs, tif.pages))
            assert len(set(array_lens)) == 1, array_lens
            for channel_idx, frame_idx, page in zip(channel_idxs, frame_idxs, tif.pages):
                summed_images[frame_idx][channel_idx] += page.asarray()

            # Add images to hdf5
            dset_name = hdf5_tools.dset_name_given_coords(Major_axis_pos, minor_axis_pos)
            for frame_idx in frame_idxs:
                h5_fpath = '{}_f{}.h5'.format(h5_bname, frame_idx)
                with h5py.File(h5_fpath, 'a') as f:
                    for idx, channel_name in enumerate(channel_names):
                        if channel_name not in f:
                            g = f.create_group(channel_name)
                        else:
                            g = f[channel_name]
        
                        out_im = np.flipud(summed_images[idx])
                        dset = g.create_dataset(dset_name, out_im.shape, dtype=out_im.dtype)
                        dset[...] = out_im
    print
            

if __name__ == '__main__':
    usg_fmt = '{} <Major_axis_idx(0,1)> <hdf5_bname> <tif_fpaths>'.format(sys.argv[0])
    if len(sys.argv) < len(usg_fmt.split()):
        sys.exit(usg_fmt)

    Major_axis_idx = int(sys.argv[1])
    assert Major_axis_idx in [0, 1], 'Major axis must be in (0, 1)'

    h5_bname = sys.argv[2]

    tif_fpaths = sys.argv[3:]

    tif_dir_to_hdf5(Major_axis_idx, h5_bname, tif_fpaths)
