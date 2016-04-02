from abc import ABCMeta, abstractmethod
import logging
#TODO: set attribute methods to be transactions
#TODO: set descriptors to forbid the word predicted
#TODO: input logging options into the settings files

logger = logging.getLogger(__name__)

class AbstractModelVisitor(object):
    __metaclass__ = ABCMeta

    maxResponseCount = None

    def __init__(self, statsModel):
        self.statsModel = statsModel

    @abstractmethod
    def train(self, verbose=False):
        """A function meant to be overridden by actual ModelVisitor classes.
              The `_train` method should prepare the machine learning model for
              classification and save that model if necessary."""

    @abstractmethod
    def predict(self, reactions, descriptorHeaders, verbose=False):
        """Return a dictionary where the key is the response descriptor being
           predicted and the value is list of prediction tuples for that
           response where the ith prediction corresponds to the ith reaction.
           EG: {<NumRxnDescriptor> "outcome" }:[(<rxn1>, 1), (<rxn2>, 2), (<rxn3>, 1), (<rxn4>, 1)]}"""
