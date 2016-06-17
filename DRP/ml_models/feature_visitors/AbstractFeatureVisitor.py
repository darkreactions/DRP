"""An abstracted version of the feature visitor."""
from abc import ABCMeta, abstractmethod
import logging

logger = logging.getLogger(__name__)


class AbstractFeatureVisitor(object):

    """The Abstracted feature visitor."""

    __metaclass__ = ABCMeta

    @abstractmethod
    def train(self, verbose=False):
        """
        A function meant to be overridden by actual FeatureVisitor classes.

        The `train` method should prepare the feature selection model for
        attribute selection and save that model if necessary.
        """
