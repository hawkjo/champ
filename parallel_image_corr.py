import os
import itertools
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import local_config
import image_processing_tools
from scipy import ndimage
from pathos.multiprocessing import ProcessingPool


project_name = 'SA15066'

def load_and_cap_stch_im(norm=True, nones=False):
    # Loading, Capping, and Shifting Image Data
    dir_60x = os.path.join(local_config.jah_04_26_dir, '60X')
    stch_60x_im = np.load(os.path.join(dir_60x, 'stitched_all.npy'))
    
    pct_995_60x = np.percentile(stch_60x_im, 99.5)
    print '60x cap:', pct_995_60x
    
    stch_60x_im[stch_60x_im > pct_995_60x] = pct_995_60x
    if not norm and not nones:
        return stch_60x_im
    
    stch_60x_im_None_pad = deepcopy(stch_60x_im)
    stch_60x_im_None_pad[stch_60x_im == 0] = None
    if not norm and nones:
        return stch_60x_im_None_pad
    
    s60norm = stch_60x_im_None_pad - 1

    if norm and nones:
        return s60norm

    s60norm[np.isnan(s60norm)] = 0
    if norm and not nones:
        return s60norm
    
    
def make_get_rotated_chunk(offset=0):
    # Extracting good chunks at various rotations
    s60norm = load_and_cap_stch_im(nones=True)
    s = 575  # Side length of chunk
    supers = 1000  # Side length of super chunk
    tl = np.where(np.isnan(s60norm[0, :]) == False)[0][0] + offset
    bl = 2200
    superchunk = s60norm[bl-supers:bl, tl:tl+supers]
    
    def get_rotated_chunk(degrees):
        rot_superchunk = ndimage.interpolation.rotate(superchunk, degrees, mode='constant', cval=float('-inf'))
        r = (rot_superchunk.shape[0] - s)/2
        c = (rot_superchunk.shape[1] - s)/2
        chunk = rot_superchunk[r:r+s, c:c+s]
        assert not np.any(chunk == float('-inf')) and chunk.shape == (s,s)
        return chunk
    return get_rotated_chunk


def make_poolable_fic_func(chunk_offset=0):
    get_rotated_chunk = make_get_rotated_chunk(chunk_offset)

    def correlate_given_rot_and_fq_w((degrees, fq_w)):
        chunk = get_rotated_chunk(degrees)
        fic = image_processing_tools.FastqImageCorrelator(project_name)
        fic.load_phiX()
        fic.set_image_data_from_ndarray(chunk, 60)
        fic.w_fq_tile = fq_w
        fic.fft_align(processors=1, verbose=False)
        return fic.max_corrs
    return correlate_given_rot_and_fq_w


def parallel_corr_calculation(degrees, fq_ws, chunk_offset=0, processors=15):
    f = make_poolable_fic_func()
    args = itertools.product(degrees, fq_ws)
    p = ProcessingPool(processors)
    results = p.map(f, args)
    return args, results


