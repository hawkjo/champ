import matplotlib
matplotlib.use('agg')
import sys
import os
import fastqimagecorrelator
import matplotlib.pyplot as plt
from pathos.multiprocessing import ProcessingPool
import local_config
import nd2reader


def get_align_params(align_param_fpath):
    #d = {name: value for line in open(align_param_fpath) for name, value in line.strip().split()}
    d = {}
    for line in open(align_param_fpath):
        if not line.strip():
            continue
        name, value = line.strip().split()
        d[name] = value

    possible_tile_keys = ['lane1tile{0}'.format(tile_num) 
                          for tile_num in map(int, d['possible_tiles'].split(','))]

    return (d['project_name'],
            d['aligning_read_names_fpath'],
            d['all_read_names_fpath'],
            int(d['objective']),
            float(d['rotation_est']),
            float(d['fq_w_est']),
            possible_tile_keys,
           )


def process_fig(align_run_name, nd2_fpath, align_param_fpath, im_idx):
    im_idx = int(im_idx)
    project_name, aligning_read_names_fpath, all_read_names_fpath, objective, rotation_est, fq_w_est, \
            possible_tile_keys = get_align_params(align_param_fpath)
    nd2 = nd2reader.Nd2(nd2_fpath)
    bname = os.path.splitext(os.path.basename(nd2_fpath))[0]
    sexcat_fpath = os.path.join(os.path.splitext(nd2_fpath)[0], '%d.cat' % im_idx)
    
    fic = fastqimagecorrelator.FastqImageCorrelator(project_name)
    fic.load_phiX()
    fic.set_image_data(im=nd2[im_idx].data, objective=objective, fpath=str(im_idx), median_normalize=True)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.align(possible_tile_keys, rotation_est, fq_w_est, hit_type=['exclusive', 'good_mutual'])
    print project_name, bname, im_idx, ','.join(tile.key for tile in fic.hitting_tiles)
    
    fig_dir = os.path.join(local_config.fig_dir, align_run_name, bname)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    ax = fic.plot_all_hits()
    ax.figure.savefig(os.path.join(fig_dir, '{}_all_hits.pdf'.format(im_idx)))

    ax = fic.plot_hit_hists()
    ax.figure.savefig(os.path.join(fig_dir, '{}_hit_hists.pdf'.format(im_idx)))

    results_dir = os.path.join(local_config.base_dir, 'results', align_run_name)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    intensity_fpath = os.path.join(results_dir, '{}_intensities.txt'.format(im_idx))
    stats_fpath = os.path.join(results_dir, '{}_stats.txt'.format(im_idx))
    fic.output_intensity_results(intensity_fpath)
    fic.write_alignment_stats(stats_fpath)

    all_fic = fastqimagecorrelator.FastqImageCorrelator(project_name)
    tile_data = local_config.fastq_tiles_given_read_name_fpath(all_read_names_fpath)
    all_fic.all_reads_fic_from_aligned_fic(fic, tile_data)
    all_read_rcs_fpath = os.path.join(results_dir, '{}_all_read_rcs.txt'.format(im_idx))
    all_fic.write_read_names_rcs(all_read_rcs_fpath)

if __name__ == '__main__':
    fmt = '{0} <align_run_name> <nd2_fpath> <align_param_file> <im_idx>'.format(sys.argv[0])
    if len(sys.argv) != len(fmt.split()):
        sys.exit('Usage: ' + fmt)
    process_fig(*sys.argv[1:])
