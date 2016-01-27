import matplotlib
matplotlib.use('agg')
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import random
import tifffile
import fastqimagecorrelator
import local_config


project_name = 'SA15128'
dname = '/home/hawkjo/BillProjects/constellation/data/20150706'
before_dir = os.path.join(dname, '20150706-SA15128_short_atto_P_and_peptide')
after_dir = os.path.join(dname, '20150706-SA15128_short_atto_P_and_peptide_NaOH')
this_fig_dir = os.path.join(local_config.fig_dir, '20150706')
this_results_dir = os.path.join(local_config.base_dir, 'results', '20150706')
for d in [this_fig_dir, this_results_dir]:
    if not os.path.isdir(d):
        os.mkdir(d)


class ConnectedComponentsFinder:
    def __init__(self, rcs, indices, conn_dist=10):
        self.rcs = rcs
        self.indices = indices
        self.conn_dist = conn_dist
        self.connect_components()

    def are_connected(self, i, j):
        return np.linalg.norm(self.rcs[self.indices[i]] - self.rcs[self.indices[j]]) <= self.conn_dist

    def connect_components(self):
        self.components = range(len(self.indices))
        for i in range(len(self.indices)):
            for j in range(i+1, len(self.indices)):
                comp_i, comp_j = self.components[i], self.components[j]
                if comp_i == comp_j:
                    continue
                elif self.are_connected(i, j):
                    self.components = [comp if comp != comp_j else comp_i for comp in self.components]

        self.connected_components = []
        for comp in set(self.components):
            self.connected_components.append(
                [idx for i, idx in enumerate(self.indices) if self.components[i] == comp]
            )


def median_normalize(im):
    med = np.median(im)
    im = im / float(med)
    im -= 1
    return im


def extract_bright_spots(im_idx):
    bname = 'xy%03dc2' % im_idx
    before_tif_fpath = os.path.join(before_dir, bname + '.tif')
    after_tif_fpath = os.path.join(after_dir, bname + '.tif')
    before_im = median_normalize(tifffile.imread(before_tif_fpath))
    after_im = median_normalize(tifffile.imread(after_tif_fpath))
    
    sexcat_fpath = os.path.join(before_dir, bname + '.cat')
    possible_tile_keys = ['lane1tile%d' % i for i in range(2106, 2109)]
    rotation_est = -0.5 - 180
    fq_w_est = 935
    
    fic = fastqimagecorrelator.FastqImageAligner(project_name)
    fic.load_phiX()
    fic.set_image_data(before_tif_fpath, 60, median_normalize=True)
    fic.set_sexcat_from_file(sexcat_fpath)
    fic.align(possible_tile_keys, rotation_est, fq_w_est=fq_w_est, snr_thresh=1.2)
    
    all_fic = fastqimagecorrelator.FastqImageAligner(project_name)
    all_fic.all_reads_fic_from_aligned_fic(fic)
    
    thresh = 1.5
    side_px = 2
    im = after_im
    offset = [-4, 1]  # manually determined
    
    all_aligned_rcs = all_fic.hitting_tiles[0].aligned_rcs
    all_aligned_rcs = all_aligned_rcs + np.tile(offset, (all_aligned_rcs.shape[0], 1))
    all_read_names = all_fic.hitting_tiles[0].read_names

    hit_indices = []
    for i, pt in enumerate(all_aligned_rcs):
        rl, ru = max(0, pt[0]-side_px), min(im.shape[0], pt[0]+side_px+1)
        cl, cu = max(0, pt[1]-side_px), min(im.shape[1], pt[1]+side_px+1)
        if 0 <= pt[0] < im.shape[0] and 0 <= pt[1] < im.shape[1] and im[rl:ru, cl:cu].max() > thresh:
            hit_indices.append(i)
    
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.matshow(im, cmap=plt.get_cmap('Blues'), vmin=0, vmax=5)
    ccf = ConnectedComponentsFinder(all_aligned_rcs, hit_indices)
    component_read_names = []
    for component in ccf.connected_components:
        color = [random.random() for _ in range(3)]
        component_read_names.append([all_read_names[idx] for idx in component])
        for idx in component:
            pt = all_aligned_rcs[idx]
            ax.plot(pt[1], pt[0], 'o', color=color, alpha=0.3)
    ax.set_title(bname)

    fig.savefig(os.path.join(this_fig_dir, bname + '.png'))
    with open(os.path.join(this_results_dir, bname + '_bright_spot_reads.txt'), 'w') as out:
        out.write('\n'.join('\t'.join(name for name in comp) for comp in component_read_names))


if __name__ == '__main__':
    usage_fmt = '%s <im_idx>' % sys.argv[0]
    if len(sys.argv) != len(usage_fmt.split()):
        sys.exit('Usage: ' + usage_fmt)
    im_idx = int(sys.argv[1])
    extract_bright_spots(im_idx)
