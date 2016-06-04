from abc import ABCMeta, abstractmethod
import DRP
import warnings


class AbstractSplitter(object):
    __metaclass__ = ABCMeta

    def __init__(self, namingStub):
        self.namingStub = namingStub
        self.namingCounter = 0

    @abstractmethod
    def split(self, data, verbose=False):
        if data.count() < 20:
            warnings.warn('You are only using {} reactions. This may cause problems (e.g. Weka SVMs require at least 10 data points in the training set)'.format(data.count()))

    def package(self, data):
        dataSet = DRP.models.DataSet.create('{}_{}'.format(self.namingStub, self.namingCounter), data)
        self.namingCounter += 1

        return dataSet
