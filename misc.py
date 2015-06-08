"""
A space for miscellaneous useful functions.
"""
import numpy as np


def next_power_of_2(x):
    return 1<<(int(np.ceil(x))-1).bit_length()


def max_2d_idx(a):
    return np.unravel_index(a.argmax(), a.shape)


def pad_to_size(M, size):
    assert len(size) == 2, 'Row and column sizes needed.'
    left_to_pad = size - np.array(M.shape) 
    return np.pad(M, ((0, left_to_pad[0]), (0, left_to_pad[1])), mode='constant')


def right_rotation_matrix(angle, degrees=True):
    if degrees:
        angle *= np.pi/180.0
    sina = np.sin(angle)
    cosa = np.cos(angle)
    return np.array([[cosa, sina],
                     [-sina, cosa]])

def rcs_given_read_names(read_names):
    return np.array([map(int, name.split(':')[-2:]) for name in read_names])
