class MicroscopeData(object):
    """
    Holds raw microscope image data and the respective Source Extractor coordinates.

    """
    def __init__(self, image, sextraction, row, column):
        self._image = image
        self._sextraction = sextraction
        self.row = row
        self.column = column

    @property
    def image(self):
        """
        Raw image from a TIRF microscope.

        """
        return self._image

    @property
    def rcs(self):
        """
        Coordinates of points as determined by Source Extractor.

        """
        return self._sextraction.rcs

    @property
    def simage(self):
        # just for fun
        return self._sextraction.image
