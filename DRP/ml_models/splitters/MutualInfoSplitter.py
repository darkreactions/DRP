from AbstractSplitter import AbstractSplitter
from DRP.models import Compound

class Splitter(AbstractSplitter):
  def split(self, reactions):

    # Create reactant-combination keys for each reactions entry.
    #reactionsKeys = create_reactant_keys(reactions)

    #TODO: implement this:
    # Partitions the reactions/keys into separate test/training reactionssets.
    #train, test = create_test_and_train_lists(reactions, reactionsKeys)

    train = reactions[:2]
    test = reactions[2:]

    return train, test


  # Returns two separate lists of data entries (a test and a training list)
  #   given data entries and the corresponding keys (in the same order).
  def create_test_and_train_lists(self, data, keys):
    test, train = [], []

    # Choose which key should go where (ie: whether a datum will
    #   be thrown to the test or training lists).
    key_in_test_map = self.build_key_in_test_map(keys)

    for i in xrange(len(keys)):
      datum = data[i]
      key = keys[i]

      if key_in_test_map[key]:
        test.append(datum)
      else:
        train.append(datum)

    return train. test

  def build_key_in_test_map(self, keys):
    key_counts = {k: keys.count(k) for k in set(keys)}
    key_list = key_counts.keys()
    key_in_test = dict()
    test_total = 0
    max_test = int(self.TEST_PERCENT*len(keys))
    for key in key_list:
      if test_total < max_test and key_counts[key] <= self.MAX_PARTITION_SIZE:
        key_in_test[key] = True
        test_total += key_counts[key]
      else:
        key_in_test[key] = False
    return key_in_test

  def count_compound_sets(self, reactions):
    compound_sets = [frozenset(Compound.objects.filter(reaction=rxn)) for rxn in reactions]
    counts = {b: sum(a==b for a in compound_sets) for b in compound_sets}
    return counts

