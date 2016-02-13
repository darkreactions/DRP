#!/usr/bin/env python
from abc import ABCMeta, abstractmethod
import logging
#TODO: set attribute methods to be transactions
#TODO: set descriptors to forbid the word predicted
#TODO: input logging options into the settings files

logger = logging.getLogger(__name__)

class AbstractModelVisitor(object):
    __metaclass__ = ABCMeta

    maxResponseCount = None
    
    distance_function = None

    def __init__(self):
        pass

    @abstractmethod
    def train(self, reactions, descriptorHeaders, filePath):
        """A function meant to be overridden by actual DistanceLeaner classes.
        The `_train` method should determine the distance function
        and save that model if necessary."""
