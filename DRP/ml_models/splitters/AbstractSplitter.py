from abc import ABCMeta, abstractmethod

class AbstractSplitter(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def split(self, data):
        pass

