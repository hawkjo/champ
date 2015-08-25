import os
import sys
import numpy as np
import nd2reader


def nd2_to_txt(nd2_fpath):
    dname = os.path.splitext(nd2_fpath)[0]
    os.mkdir(dname)

    nd2 = nd2reader.Nd2(nd2_fpath)
    for i, image in enumerate(nd2):
        out_fpath = os.path.join(dname, '%d.txt' % i)
        np.savetxt(out_fpath, image.data, fmt='%d', delimiter='\t')

if __name__ == '__main__':
    usage_fmt = '%s <nd2_file>' % sys.argv[0]
    if len(usage_fmt.split()) != len(sys.argv):
        sys.exit('Usage: ' + usage_fmt)

    nd2_fpath = sys.argv[1]
    nd2_to_txt(nd2_fpath)
