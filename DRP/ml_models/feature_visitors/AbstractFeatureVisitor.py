from abc import ABCMeta, abstractmethod
import logging

logger = logging.getLogger(__name__)

class AbstractFeatureVisitor(object):
    __metaclass__ = ABCMeta

    def __init__(self, statsModel):
        self.statsModel = statsModel

    @abstractmethod
    def train(self, reactions, descriptorHeaders, filePath):
        """A function meant to be overridden by actual FeatureVisitor classes.
              The `_train` method should prepare the feature selection model for
              attribute selection and save that model if necessary."""

    # @abstractmethod
    # def predict(self, reactions, descriptorHeaders):
    #     """Return a dictionary where the key is the response descriptor being
    #        predicted and the value is an ordered list of predictions for that
    #        response where the ith prediction corresponds to the ith reaction.
    #        EG: {<NumRxnDescriptor> "outcome" }:[1,2,1,1]}"""
