from model.tile import TileManager
from fastq import load_classified_reads
import logging


log = logging.getLogger()
log.addHandler(logging.StreamHandler())
log.setLevel(logging.DEBUG)


if __name__ == '__main__':
    read_data = load_classified_reads('phix')
    tm = TileManager(0.266666666666667, read_data)
    print(tm.tile(19).fft)
