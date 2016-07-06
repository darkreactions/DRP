"""Splitter visitor for the case when the data should not be separated."""
from .abstractSplitter import AbstractSplitter
import random
import logging
logger = logging.getLogger(__name__)


class Splitter(AbstractSplitter):
    """The splitter visitor."""

    def __init__(self, namingStub):
        """Standard splitter initialisation."""
        super(Splitter, self).__init__(namingStub)

    def split(self, reactions, verbose=False):
        """Perform the split."""
        super(Splitter, self).split(reactions, verbose=verbose)
        if verbose:
            logger.info("Training set ({}) and no test set.".format(reactions.count()))
        splits = [(self.package(reactions), self.package([]))]

        return splits
