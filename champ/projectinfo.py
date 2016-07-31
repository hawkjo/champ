import os
import h5py


def load_channels(image_directory):
    channels = set()
    for filename in os.listdir(image_directory):
        if not filename.endswith('.h5'):
            continue
        with h5py.File(os.path.join(image_directory, filename)) as h5:
            for key in h5.keys():
                channels.add(key)
    return channels
