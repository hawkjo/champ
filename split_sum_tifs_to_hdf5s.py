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
    assert len(tif_fpaths) == len(set(tif_fpaths))
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
            frame_idxs = tif.micromanager_metadata['index_map']['frame']

            # Setup summed_images
            h, w = summary['Height'], summary['Width']
            lcol, trow = 120, 0
            rcol, brow = w, h
            mrow = h/2
            mcol = w/2
            mt_row, mb_row = mrow - 20, mrow + 20
            ml_col, mr_col = mcol - 20, mcol + 20

            px_edges_given_quarter = {
                (0, 0): ((trow, mt_row), (lcol, ml_col)),
                (0, 1): ((trow, mt_row), (mr_col, rcol)),
                (1, 0): ((mb_row, brow), (lcol, ml_col)),
                (1, 1): ((mb_row, brow), (mr_col, rcol)),
            }

            side_given_quarter = {quarter: 'left' if quarter[1] == 0 else 'right'
                                  for quarter in px_edges_given_quarter.keys()}
            quarters_given_side = {side: [] for side in side_given_quarter.values()}
            for quarter, side in side_given_quarter.items():
                quarters_given_side[side].append(quarter)

            summed_images = {frame_idx:
                             {channel_name + '_' + side:
                              {
                                  quarter: np.zeros((t-b, r-l), dtype=np.int)
                                  for quarter, ((b, t), (l, r)) in px_edges_given_quarter.items()
                                  if quarter in side_quarters
                              }
                              for channel_name in set(channel_names)
                              for side, side_quarters in quarters_given_side.items()
                             }
                             for frame_idx in set(frame_idxs)
                            }

            # Add images
            array_lens = map(len, (channel_idxs, frame_idxs, tif.pages))
            assert len(set(array_lens)) == 1, array_lens
            for channel_idx, frame_idx, page in zip(channel_idxs, frame_idxs, tif.pages):
                for quarter, ((b, t), (l, r)) in px_edges_given_quarter.items():
                    sided_channel_name = channel_names[channel_idx] + '_' + side_given_quarter[quarter]
                    summed_images[frame_idx][sided_channel_name][quarter] += page.asarray()[b:t, l:r]

            # Add images to hdf5
            for frame_idx in set(frame_idxs):
                h5_fpath = '{}_f{}.h5'.format(h5_bname, frame_idx)
                with h5py.File(h5_fpath, 'a') as f:
                    for channel_idx, channel_name in enumerate(channel_names):
                        for quarter in px_edges_given_quarter.keys():
                            sided_channel_name = channel_name + '_' + side_given_quarter[quarter]
                            if sided_channel_name not in f:
                                g = f.create_group(sided_channel_name)
                            else:
                                g = f[sided_channel_name]

                            qr, qc = quarter
                            loc_Major_axis_pos = 2*Major_axis_pos + (1-qr)
                            loc_minor_axis_pos = 2*minor_axis_pos + (1-qc)
                            dset_name = hdf5_tools.dset_name_given_coords(loc_Major_axis_pos,
                                                                          loc_minor_axis_pos)
                            out_im = np.flipud(summed_images[frame_idx][sided_channel_name][quarter])
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
