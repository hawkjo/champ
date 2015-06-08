import os
import itertools
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import local_config
import fastqimagecorrelator
from scipy import ndimage
from pathos.multiprocessing import ProcessingPool


def load_and_cap_stch_im(project_name, norm=True, nones=False):
    # Loading, Capping, and Shifting Image Data
    if project_name == 'SA15066':
        im_dir = os.path.join(local_config.jah_04_26_dir, '60X')
    elif project_name == 'SA15097':
        im_dir = os.path.join(local_config.data_dir, '20150528/stitched_images/SA15097_before')
    else:
        raise ValueError('Invalid project name: {0}'.format(project_name))
    stch_im = np.load(os.path.join(im_dir, 'stitched_all.npy'))
    
    pct_995 = np.percentile(stch_im, 99.5)
    print 'Cap:', pct_995
    
    stch_im[stch_im > pct_995] = pct_995
    if not norm and not nones:
        return stch_im
    
    stch_im_None_pad = deepcopy(stch_im)
    stch_im_None_pad[stch_im == 0] = None
    if not norm and nones:
        return stch_im_None_pad
    
    stch_im_norm = stch_im_None_pad - 1

    if norm and nones:
        return stch_im_norm

    stch_im_norm[np.isnan(stch_im_norm)] = 0
    if norm and not nones:
        return stch_im_norm


def make_get_rotated_chunk(project_name, offset=0):
    # Extracting good chunks at various rotations
    stch_im_norm = load_and_cap_stch_im(project_name, nones=True)
    s = 575  # Side length of chunk
    supers = 1000  # Side length of super chunk
    #tl = np.where(np.isnan(stch_im_norm[0, :]) == False)[0][0] + offset
    #bl = 2200
    tl = 200 + offset
    bl = stch_im_norm.shape[0] - 200
    superchunk = stch_im_norm[bl-supers:bl, tl:tl+supers]
    
    def get_rotated_chunk(degrees):
        rot_superchunk = ndimage.interpolation.rotate(superchunk, degrees, mode='constant', cval=float('-inf'))
        r = (rot_superchunk.shape[0] - s)/2
        c = (rot_superchunk.shape[1] - s)/2
        chunk = rot_superchunk[r:r+s, c:c+s]
        assert not np.any(chunk == float('-inf')) and chunk.shape == (s,s)
        return chunk
    return get_rotated_chunk


def make_poolable_fic_func(project_name, chunk_offset=0):
    get_rotated_chunk = make_get_rotated_chunk(project_name, chunk_offset)

    def correlate_given_rot_and_fq_w((degrees, fq_w)):
        chunk = get_rotated_chunk(degrees)
        fic = fastqimagecorrelator.FastqImageCorrelator(project_name)
        fic.load_phiX()
        fic.set_image_data(im=chunk, objective=60)
        fic.w_fq_tile = fq_w
        fic.fft_align(processors=1, verbose=False)
        return fic.max_corrs
    return correlate_given_rot_and_fq_w


def parallel_corr_calculation(project_name, degrees, fq_ws, chunk_offset=0, processors=15):
    f = make_poolable_fic_func(project_name, chunk_offset)
    args = itertools.product(degrees, fq_ws)
    p = ProcessingPool(processors)
    results = p.map(f, args)
    return args, results


