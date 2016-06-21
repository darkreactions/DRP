"""Contains an abstracted splitter class on which other splitters may be based."""
from abc import ABCMeta, abstractmethod
import DRP
import warnings


class AbstractSplitter(object):

    """An abstracted splitter class on which other splitters may be based."""

    __metaclass__ = ABCMeta

    def __init__(self, namingStub):
        """Allow us to name the datasets with a namingstub."""
        self.namingStub = namingStub
        self.namingCounter = 0

    @abstractmethod
    def split(self, data, verbose=False):
        """Actually perform the split."""
        if data.count() < 20:  # TODO: This should not be done here.
            warnings.warn(
                'You are only using {} reactions. This may cause problems (e.g. Weka SVMs require at least 10 data points in the training set)'.format(data.count()))

    def package(self, data):
        """Save the splits as datasets in the database."""
        dataSet = DRP.models.DataSet.create('{}_{}'.format(
            self.namingStub, self.namingCounter), data)
        self.namingCounter += 1

        return dataSet
