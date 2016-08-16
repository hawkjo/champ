import matplotlib.pyplot as plt
import h5py


def imshow(filename, channel, row, column):
    with h5py.File(filename, 'r') as h5:
        image = h5[channel]['(Major, minor) = ({}, {})'.format(column, row)].value
        fig = plt.figure(figsize=(15, 15))
        plt.imshow(image, cmap='viridis')
        plt.show()
