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
        dsrs = [DRP.models.DataSetRelation(dataSet=dataSet, reaction=datum) for datum in data]
        DRP.models.DataSetRelation.objects.bulk_create(dsrs)

        return dataSet
