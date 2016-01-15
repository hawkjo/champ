import os


class SEConfig(object):
    def __enter__(self):
        self._create_config_files()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._delete_config_files()

    def _delete_config_files(self):
        for filename in ('default.sex', 'spot.param', 'default.conv'):
            try:
                os.unlink(filename)
            except OSError:
                pass

    def _create_config_files(self):
        default_text = """DETECT_THRESH 2
DEBLEND_NTHRESH 64
DEBLEND_MINCONT 0.00005
"""
        with open('default.sex', 'w+') as f:
            f.write(default_text)

        spot_text = """X_IMAGE
Y_IMAGE
FLUX_AUTO
FLUXERR_AUTO
FLAGS
A_IMAGE
B_IMAGE
THETA_IMAGE
"""
        with open('spot.param', 'w+') as f:
            f.write(spot_text)

        convolution_text = """CONV NORM
1 2 1
2 4 2
1 2 1
"""
        with open('default.conv', 'w+') as f:
            f.write(convolution_text)
