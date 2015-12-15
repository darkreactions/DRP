from abc import ABCMeta, abstractmethod
import DRP

class AbstractSplitter(object):
    __metaclass__ = ABCMeta


    def __init__(self, namingStub):
        self.namingStub = namingStub
        self.namingCounter = 0

    @abstractmethod
    def split(self, data):
        pass

    def package(self, data):
        dataSet = DRP.models.DataSet(name=self.namingStub + '_{}'.format(self.namingCounter))
        dataSet.save()
        self.namingCounter+=1
        for datum in data:
            dsr = DRP.models.DataSetRelation(dataSet = dataSet, reaction=datum)
            dsr.save()
        return dataSet
