from AbstractSplitter import AbstractSplitter
import random

class Splitter(AbstractSplitter):
  def __init__(self, namingStub, test_percent=0.33):
    super(Splitter, self).__init__(namingStub)
    self.test_percent = test_percent

  def split(self, reactions, verbose=False):
    num_rxns = reactions.count()
    test_size = int(self.test_percent*num_rxns)

    # Split the reactions' IDs into two randomly-organized buckets.
    rxn_ids = [reaction.id for reaction in reactions]
    random.shuffle(rxn_ids)

    test = reactions.filter(id__in=rxn_ids[:test_size])
    train = reactions.filter(id__in=rxn_ids[test_size:])

    if verbose:
        print "Split into train ({}), test ({})".format(train.count(), test.count())
    
    splits = [ (self.package(train), self.package(test)) ]

    return splits
