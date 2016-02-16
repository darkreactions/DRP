from abc import ABCMeta, abstractmethod
import logging

logger = logging.getLogger(__name__)

class AbstractDistanceLearner(object):
    __metaclass__ = ABCMeta

    maxResponseCount = None
    
    distance_function = None

    def __init__(self):
        self.distance_function = None

    @abstractmethod
    def train(self, reactions, predictorHeaders, responseHeaders):
        """A function meant to be overridden by actual DistanceLeaner classes.
        The `_train` method should determine the distance function
        and save that model if necessary."""

    def distance_function_from_matrix(self, matrix):
        def _distance_function(x, y):
            return x.dot(matrix.dot(y))
        return _distance_function
