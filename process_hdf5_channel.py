import matplotlib
matplotlib.use('agg')
import sys
import os
import fastqimagecorrelator
import matplotlib.pyplot as plt
from pathos.multiprocessing import ProcessingPool
import local_config
import h5py
import hdf5_tools
import random
import numpy as np


def get_align_params(align_param_fpath):
    #d = {name: value for line in open(align_param_fpath) for name, value in line.strip().split()}
    d = {}
    for line in open(align_param_fpath):
        if not line.strip():
            continue
        name, value = line.strip().split()
        d[name] = value

    try:
        strategy = d['strategy']
    except:
        strategy = 'slow'
    assert strategy in ['fast', 'slow'], strategy

    try:
        snr_thresh = float(d['snr_thresh'])
    except:
        snr_thresh = 1.2

    try:
        min_hits = int(d['min_hits'])
    except:
        min_hits = 15

    return (d['project_name'],
            d['aligning_read_names_fpath'],
            d['all_read_names_fpath'],
            int(d['objective']),
            float(d['rotation_est']),
            float(d['fq_w_est']),
            int(d['min_tile_num']),
            int(d['max_tile_num']),
            strategy,
            snr_thresh,
            min_hits
           )


def fast_possible_tile_keys(Maj_pos, min_pos, max_pos, min_tile, max_tile):
    expected_tile = int(min_tile + float(Maj_pos - min_pos)/(max_pos - min_pos) * (max_tile - min_tile))
    return tile_keys_given_nums(range(expected_tile - 1, expected_tile + 2))


def tile_keys_given_nums(tile_nums):
    return ['lane1tile{0}'.format(tile_num) for tile_num in tile_nums]


def process_hdf5(align_run_name, hdf5_fpath, align_param_fpath, channel, im_idx):
    print 'Loading basic settings...'
    im_idx = int(im_idx)
    hdf5_fname = os.path.basename(hdf5_fpath)
    hdf5_bname = os.path.splitext(hdf5_fname)[0]
    sexcat_dir = hdf5_bname + '_sex_files'
    assert os.path.isdir(sexcat_dir), 'Sextractor directory does not exist: {}'.format(sexcat_dir)

    fig_dir = os.path.join(local_config.fig_dir, align_run_name, hdf5_bname)
    results_dir = os.path.join(local_config.base_dir, 'results', align_run_name, hdf5_bname)
    for d in [fig_dir, results_dir]:
        if not os.path.exists(d):
            os.makedirs(d)

    project_name, aligning_read_names_fpath, all_read_names_fpath, objective, rotation_est, fq_w_est, \
            min_tile_num, max_tile_num, strategy, snr_thresh, min_hits = get_align_params(align_param_fpath)

    with h5py.File(hdf5_fpath) as f:
        g = f[channel]

        print 'Setting possible_tile_keys_func...'
        all_Major_pos = map(hdf5_tools.Major_pos_given_dset_name, g.keys())
        max_pos, min_pos = max(all_Major_pos), min(all_Major_pos)
        if strategy == 'fast':
            possible_tile_keys_func = lambda x: fast_possible_tile_keys(x, min_pos, max_pos, min_tile_num, max_tile_num)
        elif strategy == 'slow':
            possible_tile_keys_func = lambda _: tile_keys_given_nums(range(min_tile_num, max_tile_num+1))

        def process_im_wrapper(tup):
            dset_name, im = tup
            im_bname = hdf5_tools.bname_given_channel_and_dset_name(channel, dset_name)
            Major_pos = hdf5_tools.Major_pos_given_dset_name(dset_name)
            possible_tile_keys = possible_tile_keys_func(Major_pos)
            process_im(im,
                       im_bname,
                       sexcat_dir,
                       fig_dir,
                       results_dir,
                       possible_tile_keys,
                       align_param_fpath)

        dset_names = g.keys()
        dset_names.sort()
        dset_name = dset_names[im_idx]
        im = np.array(g[dset_name])
        process_im_wrapper((dset_name, im))

        return 
#        input_params = [(dset_name, np.array(g[dset_name]), tile_data, all_tile_data) for dset_name in g.keys()]
#        print 'Processing images...'
#        p = ProcessingPool(num_threads)
#        p.map(process_im_wrapper, input_params[:20]).get(1234567)
#        p.close()

        
def process_im(im,
               im_bname,
               sexcat_dir,
               fig_dir,
               results_dir,
               possible_tile_keys,
               align_param_fpath):
    sexcat_fpath = os.path.join(sexcat_dir, '{}.cat'.format(im_bname))
    
    all_read_rcs_fpath = os.path.join(results_dir, '{}_all_read_rcs.txt'.format(im_bname))

    if os.path.isfile(all_read_rcs_fpath):
        print im_bname, 'already done.'

    project_name, aligning_read_names_fpath, all_read_names_fpath, objective, rotation_est, fq_w_est, \
            min_tile_num, max_tile_num, strategy, snr_thresh, min_hits = get_align_params(align_param_fpath)

    intensity_fpath = os.path.join(results_dir, '{}_intensities.txt'.format(im_bname))
    stats_fpath = os.path.join(results_dir, '{}_stats.txt'.format(im_bname))
    fic = fastqimagecorrelator.FastqImageCorrelator(project_name)
    tile_data=local_config.fastq_tiles_given_read_name_fpath(aligning_read_names_fpath)
    fic.load_reads(tile_data)
    fic.set_image_data(im=im, objective=objective, fpath=im_bname, median_normalize=True)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.align(possible_tile_keys, rotation_est, fq_w_est, snr_thresh=snr_thresh, min_hits=min_hits, hit_type=['exclusive', 'good_mutual'])
    print project_name, im_bname, ','.join(tile.key for tile in fic.hitting_tiles)
    
    fic.output_intensity_results(intensity_fpath)
    fic.write_alignment_stats(stats_fpath)

    ax = fic.plot_all_hits()
    ax.figure.savefig(os.path.join(fig_dir, '{}_all_hits.pdf'.format(im_bname)))

    ax = fic.plot_hit_hists()
    ax.figure.savefig(os.path.join(fig_dir, '{}_hit_hists.pdf'.format(im_bname)))

    all_fic = fastqimagecorrelator.FastqImageCorrelator(project_name)
    all_tile_data = local_config.fastq_tiles_given_read_name_fpath(all_read_names_fpath)
    all_fic.all_reads_fic_from_aligned_fic(fic, all_tile_data)
    all_fic.write_read_names_rcs(all_read_rcs_fpath)

    print
    print '!'*80
    print 
    print 'Alignment found:', im_bname
    print
    print '!'*80
    print

if __name__ == '__main__':
    fmt = '{0} <align_run_name> <hdf5_fpath> <align_param_file> <channel> <im_idx>'.format(sys.argv[0])
    if len(sys.argv) != len(fmt.split()):
        sys.exit('Usage: ' + fmt)
    process_hdf5(*sys.argv[1:])
