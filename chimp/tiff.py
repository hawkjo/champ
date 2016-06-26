from abc import abstractproperty, abstractmethod


class BaseTifStack(object):
    @abstractproperty
    def axes(self):
        raise NotImplementedError

    @abstractmethod
    def __iter__(self):
        raise NotImplementedError


class TifsPerFieldOfView(BaseTifStack):
    pass


class TifsPerConcentration(BaseTifStack):
    pass


class TIFSingleFieldOfView(object):
    pass