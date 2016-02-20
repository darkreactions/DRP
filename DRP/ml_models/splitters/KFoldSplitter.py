from AbstractSplitter import AbstractSplitter
import random

class Splitter(AbstractSplitter):
  def __init__(self, namingStub):
    super(Splitter, self).__init__(namingStub)
    self.k = 4

  def split(self, reactions):
    num_rxns = reactions.count()
    bucket_size = num_rxns/self.k

    # Split the reactions' IDs into K randomly-organized buckets.
    rxn_ids = [reaction.id for reaction in reactions]
    random.shuffle(rxn_ids)
    buckets = [rxn_ids[i:i+bucket_size] for i in range(0, len(rxn_ids), bucket_size)]

    # Create a Reaction QuerySet for the training and test buckets every fold.
    splits = []
    for i in range(self.k):
      train = reactions.filter(id__in=[item for b in buckets[:i]+buckets[i+1:]
                                            for item in b])
      test = reactions.filter(id__in=buckets[i])
      splits.append( (self.package(train), self.package(test)) )

    return splits
