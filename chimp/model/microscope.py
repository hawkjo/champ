class MicroscopeData(object):
    """
    Holds raw microscope image data and the respective Source Extractor coordinates.

    """
    def __init__(self, image, sextraction, row, column):
        self._image = image
        self.column = column
        self.row = row
        self.shape = image.shape
        self._sextraction = sextraction

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
