class GridImages(object):
    # TODO: This should return MicroscopeData instead of raw ND2 Images

    def __init__(self, nd2, alignment_channel=None, channel_offset=None):
        """
        Since some ND2s were created where multiple channels had the same name, we can't always use the channel name, though we will be able to going forward
        now that we saw that that was happening.

        """
        self._nd2 = nd2
        self._rows = None
        self._columns = None
        if alignment_channel is not None:
            self._channel_offset = self._determine_channel_offset(alignment_channel)
        elif channel_offset is not None:
            self._channel_offset = channel_offset
        else:
            raise ValueError("You need to provide a name or an integer offset for the alignment channel.")
        self._parse_grid()

    def left_iter(self):

        for column in range(self._width):
            for row in range(self._height):
                index = self._get_indexes(row, column)[self._channel_offset]
                yield self._nd2[index]

    def _determine_channel_offset(self, channel_name):
        maximum_reasonable_number_of_channels = 10
        for n, image in zip(range(maximum_reasonable_number_of_channels), self._nd2):
            if image.channel == channel_name:
                return n
        raise ValueError('The channel you set to be used for alignment was not found in the given ND2.')

    def _get_first_offset_number(self, row, column):
        return column * len(self._nd2.channels) + (row * self._width)

    def _get_indexes(self, row, column):
        first = self._get_first_offset_number(row, column)
        return tuple((first + i for i in range(len(self._nd2.channels))))

    def _parse_grid(self):
        try:
            coord_info = self._nd2._parser.raw_metadata.image_metadata['SLxExperiment']['uLoopPars']['Points']['']
        except (KeyError, IndexError):
            coord_info = self._nd2._parser.raw_metadata.image_metadata['SLxExperiment']['ppNextLevelEx']['']['uLoopPars']['Points']['']
        pos_names = [pt['dPosName'] for pt in coord_info]
        self._height = len(set(name[0] for name in pos_names))
        self._width = len(set(name[1:] for name in pos_names)) * len(self._nd2.channels)
