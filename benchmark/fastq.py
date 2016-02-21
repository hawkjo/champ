from chimp.fastq import load_tiles
import time


def measure(f):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = f(*args, **kwargs)
        print("%s took %s seconds" % (f.__name__, round(time.time() - start, 2)))
        return result
    return wrapper


@measure
def load_fastq_data():
    tiles = list(load_tiles("/home/jim/Desktop/ngs/SA15243/all_fastqs"))
    print(tiles)
    print(len(tiles))

if __name__ == '__main__':
    load_fastq_data()
