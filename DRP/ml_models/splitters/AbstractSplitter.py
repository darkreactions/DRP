from abc import ABCMeta, abstractmethod
import DRP

class AbstractSplitter(object):
    __metaclass__ = ABCMeta

    def __init__(self, namingStub):
        self.namingStub = namingStub
        self.namingCounter = 0

    @abstractmethod
    def split(self, data, verbose=False):
        pass

    def package(self, data):
        dataSet = DRP.models.DataSet.create('{}_{}'.format(self.namingStub, self.namingCounter), data)
        self.namingCounter += 1
    
        return dataSet
