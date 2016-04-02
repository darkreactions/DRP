from AbstractSplitter import AbstractSplitter
import random

class Splitter(AbstractSplitter):
    def __init__(self, namingStub, num_folds=4):
        super(Splitter, self).__init__(namingStub)
        self.k = num_folds

    def split(self, reactions, verbose=False):
        # Split the reactions' IDs into K randomly-organized buckets.
        rxn_ids = [reaction.id for reaction in reactions]
        random.shuffle(rxn_ids)
        buckets = [rxn_ids[i::self.k] for i in xrange(self.k)]

        if verbose:
            print "Split into {} buckets with sizes: {}".format(len(buckets), map(len, buckets))
            
        splits = []
        for i in range(self.k):
            train = reactions.filter(id__in=[item for b in buckets[:i]+buckets[i+1:]
                                        for item in b])
            test = reactions.filter(id__in=buckets[i])
            splits.append( (self.package(train), self.package(test)) )


        return splits
