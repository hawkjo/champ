import numpy as np


def next_power_of_2(x):
    return 1 << (int(np.ceil(x)) - 1).bit_length()


def calculate_pad_size(tile_height, tile_width, image_height, image_width):
    return next_power_of_2(max(tile_height, image_height, tile_width, image_width))


def pad_image(image, pad_to_size):
    pad = pad_to_size - image.shape[0], pad_to_size - image.shape[1]
    return np.pad(image, ((0, pad[0]), (0, pad[1])), mode='constant')
