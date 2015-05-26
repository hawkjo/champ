import os
import fastqimagecorrelator
import matplotlib.pyplot as plt
from pathos.multiprocessing import ProcessingPool
import local_config


project_name = 'SA15066'
objective = 40

def process_fig(i):
    bname = 'xy%03dc1' % i
    print bname
    data_dir = os.path.join(local_config.base_dir, 'data/2015_04_26_imaging_data/%dX' % objective)
    tif_fpath = os.path.join(data_dir, 'raw', bname + '.tif')
    sexcat_fpath = os.path.join(data_dir, 'sextractor_raw', bname + '.cat')
    
    fic = fastqimagecorrelator.FastqImageCorrelator(project_name)
    if objective == 40:
        fic.w_fq_tile = 917
    fic.load_phiX()
    fic.set_image_data(tif_fpath, objective, median_normalize=True)
    fic.set_fastq_tile_mappings()
    fic.set_all_fastq_image_data()
    fic.rotate_all_fastq_data(-2.25-180)
    fic.find_hitting_tiles(['lane1tile%d' % i for i in range(2106, 2114)])
    print [tile.key for tile in fic.hitting_tiles]
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.find_hits(good_mutual_hit_thresh=5)
    
    fic.plot_all_hits()
    plt.savefig(os.path.join(local_config.fig_dir,
                             '04_26',
                             '%dX' % objective,
                             '%s_all_hits.pdf' % bname))
    plt.close()
    fic.plot_hit_hists()
    plt.savefig(os.path.join(local_config.fig_dir,
                             '04_26',
                             '%dX' % objective,
                             '%s_hit_hists.pdf' % bname))
    plt.close()
    fic.plot_threshold_gmm(force=True)
    plt.savefig(os.path.join(local_config.fig_dir,
                             '04_26',
                             '%dX' % objective,
                             '%s_gmm_thresh.pdf' % bname))
    plt.close()
    
    name_list_fpath = os.path.join(local_config.base_dir,
            'data',
            '04_26_valid_names_%dX.txt' % objective)
    out_fpath = os.path.join(local_config.base_dir,
                             'results',
                             '04_26',
                             '%dX' % objective,
                             '%s_intensities.txt' % bname)
    fic.extract_intensity_from_name_list(name_list_fpath, out_fpath)


if __name__ == '__main__':
    import sys
    i = int(sys.argv[1])
    process_fig(i)
