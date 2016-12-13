import matplotlib.pyplot as plt
import h5py
from champ.clusters import Clusters
import os
from champ.grid import GridImages
from champ import align
from champ.fastqimagealigner import FastqImageAligner
from champ.chip import Miseq


def cluster_filepath(filename, image):
    return os.path.join(basename(filename), '%s.cat' % image.index)


def basename(filename):
    return os.path.splitext(filename)[0]


def load_image(filename, channel, row, column):
    with h5py.File(filename, 'r') as h5:
        grid = GridImages(h5, channel)
        return grid.get(row, column)


def imshow(filename, channel, row, column):
    image = load_image(filename, channel, row, column)
    cluster_path = cluster_filepath(filename, image)
    with open(cluster_path) as f:
        clusters = Clusters(f)
        fig = plt.figure(figsize=(15, 15))
        plt.imshow(image, cmap='viridis')
        plt.scatter(clusters.cs(), clusters.rs(), facecolors='none', edgecolors='red', alpha=0.7)
        plt.show()


def get_fastq_image_aligner(read_names_path):
    reads = align.load_read_names(read_names_path)
    fia = FastqImageAligner(0.2666666)
    fia.load_reads(reads)
    return fia


def align_one_image(filename, channel, row, column, read_names_path):
    fia = get_fastq_image_aligner(read_names_path)
    image = load_image(filename, channel, row, column)
    cluster_path = cluster_filepath(filename, image)
    miseq_chip = Miseq(False)
    for tiles in (miseq_chip.left_side_tiles, miseq_chip.right_side_tiles):
        fia.set_image_data(image, 0.26666666666)
        fia.set_clusters_from_file(cluster_path)
        fia.rough_align(tiles, miseq_chip.rotation_estimate, miseq_chip.tile_width, snr_thresh=1.2)
        if fia.hitting_tiles:
            print("Image aligned to %s", fia.hitting_tiles)
            return fia
    print("Image did not align to any tiles!")
