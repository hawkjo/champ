class MicroscopeData(object):
    """
    Holds raw microscope image data and the respective Source Extractor coordinates.

    """
    def __init__(self, image, sextraction):
        self._image = image
        self._sextraction = sextraction

    @property
    def image(self):
        return self._image

    @property
    def rcs(self):
        return self._sextraction.rcs
