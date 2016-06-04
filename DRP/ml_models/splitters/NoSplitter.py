from AbstractSplitter import AbstractSplitter
import random


class Splitter(AbstractSplitter):

    def __init__(self, namingStub):
        super(Splitter, self).__init__(namingStub)

    def split(self, reactions, verbose=False):
        super(Splitter, self).split(reactions, verbose=verbose)
        if verbose:
            print "Training set ({}) and no test set.".format(reactions.count())
        splits = [(self.package(reactions), self.package([]))]

        return splits
