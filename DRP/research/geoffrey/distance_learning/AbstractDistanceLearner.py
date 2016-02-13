from abc import ABCMeta, abstractmethod
import logging

logger = logging.getLogger(__name__)

class AbstractDistanceLearner(object):
    __metaclass__ = ABCMeta

    maxResponseCount = None
    
    distance_function = None

    def __init__(self):
        pass

    @abstractmethod
    def train(self, reactions, predictorHeaders, responseHeaders):
        """A function meant to be overridden by actual DistanceLeaner classes.
        The `_train` method should determine the distance function
        and save that model if necessary."""
