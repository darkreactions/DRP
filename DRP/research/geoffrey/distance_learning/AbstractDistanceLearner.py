from abc import ABCMeta, abstractmethod
import logging
from cPickle import dump

logger = logging.getLogger(__name__)


class AbstractDistanceLearner(object):
    __metaclass__ = ABCMeta

    maxResponseCount = None

    @abstractmethod
    def train(self):
        """A function meant to be overridden by actual DistanceLeaner classes.
        The `train` method should determine the distance function
        and save that model if necessary."""

    @abstractmethod
    def dist(self, x, y):
        """Returns the distance between x and y under the metric"""
