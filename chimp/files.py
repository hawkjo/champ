import os


class ImageFiles(object):
    def __init__(self, filenames):
        self._filenames = filenames

    def __len__(self):
        return len(self._filenames)

    @property
    def filenames(self):
        for f in self._filenames:
            yield f

    @property
    def directories(self):
        for f in self._filenames:
            yield os.path.splitext(f)[0]


def load_image_files():
    filenames = [f for f in os.listdir(os.getcwd()) if f.endswith('.h5')]
    return ImageFiles(filenames)


def ensure_image_data_directory_exists(h5_filename):
    """
    Creates a directory based on the HDF5 filenames in order to store data derived from them.

    """
    new_directory = os.path.join(h5_filename)
    if not os.path.isdir(new_directory):
        os.mkdir(new_directory)
