from abc import ABCMeta, abstractmethod
import logging
from cPickle import dump

logger = logging.getLogger(__name__)

class AbstractDistanceLearner(object):
    __metaclass__ = ABCMeta

    maxResponseCount = None

    def __init__(self, reactions, predictorHeaders, responseHeaders):
        self.reactions = reactions
        self.predictorHeaders = predictorHeaders
        self.responseHeaders = responseHeaders

    @abstractmethod
    def train(self):
        """A function meant to be overridden by actual DistanceLeaner classes.
        The `_train` method should determine the distance function
        and save that model if necessary."""

    @abstractmethod
    def dist(self, x, y):
        """Returns the distance between x and y under the metric"""
            
    def save(self, writeable):
        dump(self, writeable)
