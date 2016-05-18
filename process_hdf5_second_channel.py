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
        min_hits = int(d['min_hits'])
    except:
        min_hits = 15

    return (d['project_name'],
            d['aligning_read_names_fpath'],
            d['all_read_names_fpath'],
            int(d['objective']),
            d['aligned_channel'],
            min_hits
           )


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

    project_name, aligning_read_names_fpath, all_read_names_fpath, objective, \
            aligned_channel, min_hits = get_align_params(align_param_fpath)

    with h5py.File(hdf5_fpath) as f:
        if channel not in f or aligned_channel not in f:
            sys.exit('Invalid channel in {}. Options: {}'.format((channel, aligned_channel),
                                                                 ', '.join(f.keys())))
        g = f[channel]
        def process_im_wrapper(tup):
            dset_name, im, im_idx = tup
            im_bname = hdf5_tools.bname_given_channel_and_dset_name(channel, dset_name)
            aligned_im_bname = hdf5_tools.bname_given_channel_and_dset_name(aligned_channel, dset_name)
            process_im(im,
                       im_idx,
                       im_bname,
                       aligned_im_bname,
                       sexcat_dir,
                       fig_dir,
                       results_dir,
                       align_param_fpath)

        dset_names = g.keys()
        dset_names.sort()
        dset_name = dset_names[im_idx]
        im = np.array(g[dset_name])
        process_im_wrapper((dset_name, im, im_idx))

        
def process_im(im,
               im_idx,
               im_bname,
               aligned_im_bname,
               sexcat_dir,
               fig_dir,
               results_dir,
               align_param_fpath):
    sexcat_fpath = os.path.join(sexcat_dir, '{}.cat'.format(im_bname))
    
    all_read_rcs_fpath = os.path.join(results_dir, '{}_all_read_rcs.txt'.format(im_bname))
    intensity_fpath = os.path.join(results_dir, '{}_intensities.txt'.format(im_bname))
    stats_fpath = os.path.join(results_dir, '{}_stats.txt'.format(im_bname))
    aligned_stats_fpath = os.path.join(results_dir, '{}_stats.txt'.format(aligned_im_bname))

    if os.path.isfile(all_read_rcs_fpath):
        print im_bname, 'already done.'
    if not os.path.isfile(aligned_stats_fpath):
        sys.exit('{}: No previous alignment.'.format(im_bname))

    project_name, aligning_read_names_fpath, all_read_names_fpath, objective, \
            aligned_channel, min_hits = get_align_params(align_param_fpath)

    fic = fastqimagecorrelator.FastqImageCorrelator(project_name)
    tile_data=local_config.fastq_tiles_given_read_name_fpath(aligning_read_names_fpath)
    fic.load_reads(tile_data)
    fic.set_image_data(im=im, objective=objective, fpath=im_bname, median_normalize=True)
    fic.set_sexcat_from_file(sexcat_fpath)

    fic.alignment_from_alignment_file(aligned_stats_fpath)
    fic.precision_align_only(min_hits=min_hits)
    print project_name, im_bname, ','.join(tile.key for tile in fic.hitting_tiles)
    
    fic.output_intensity_results(intensity_fpath)
    fic.write_alignment_stats(stats_fpath)

    if im_idx % 20 == 0:
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
