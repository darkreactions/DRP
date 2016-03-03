from abc import ABCMeta, abstractmethod
import logging

logger = logging.getLogger(__name__)

class AbstractFeatureVisitor(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def train(self, reactions, descriptorHeaders, filePath, verbose=False):
        """A function meant to be overridden by actual FeatureVisitor classes.
              The `train` method should prepare the feature selection model for
              attribute selection and save that model if necessary."""

