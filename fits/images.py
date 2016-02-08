import os


def get_nd2_filenames():
    """
    Finds the ND2s in the current directory and returns their names without the .nd2 extension.

    """
    for f in os.listdir(os.getcwd()):
        if f.endswith(".nd2"):
            yield os.path.splitext(f)[0]


def make_image_data_directory(nd2_filename):
    """
    Creates a directory based on the ND2 filenames in order to store data derived from them.

    """
    if not os.path.isdir(nd2_filename):
        os.mkdir(nd2_filename)
