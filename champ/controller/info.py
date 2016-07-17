import logging
import os
import h5py

log = logging.getLogger(__name__)


def main(clargs):
    channels = set()
    for filename in os.listdir(clargs.image_directory):
        if not filename.endswith('.h5'):
            continue
        with h5py.File(os.path.join(clargs.image_directory, filename)) as h5:
            for key in h5.keys():
                channels.add(key)
    print("Channels:")
    for channel in channels:
        print(channel)
