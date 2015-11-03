from AbstractSplitter import AbstractSplitter
from DRP.models import Compound
import random

class Splitter(AbstractSplitter):
  def __init__(self):
    super(Splitter, self).__init__()
    self.TEST_PERCENT = 0.30
    self.MAX_PARTITION_SIZE = 35 # Magic #TODO: Please de-magic this.
    self.MIN_TRAIN_SIZE = 10 # WEKA requires at least 10 training points for SVM.

  def split(self, reactions):

    # Partition the reactions based on what compounds they contain.
    key_counts = self._count_compound_sets(reactions).items()
    random.shuffle(key_counts)

    max_test_size = int(self.TEST_PERCENT*len(key_counts))
    test_size = 0
    train_size = 0

    # Determine which partitions should be tested.
    test_keys = []
    for key, count in key_counts:
      if train_size < self.MIN_TRAIN_SIZE:
        train_size += count
      elif test_size < max_test_size and count <= self.MAX_PARTITION_SIZE:
        test_keys.append(key)
        test_size += count

    # Get the primary keys of the reactions that should be tested.
    test_reactions = set()
    for key in test_keys:
      for reaction in reactions.filter(compounds__in=key):
        test_reactions.add( reaction.pk )

    # Actually perform the query.
    train = reactions.exclude(pk__in=test_reactions)
    test = reactions.filter(pk__in=test_reactions)

    return train, test

  def _count_compound_sets(self, reactions):
    compound_sets = [frozenset(Compound.objects.filter(reaction=rxn)) for rxn in reactions]
    counts = {b: sum(a==b for a in compound_sets) for b in compound_sets}
    return counts

