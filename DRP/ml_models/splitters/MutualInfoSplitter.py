from AbstractSplitter import AbstractSplitter
from DRP.models.Compound import Compound
from django.db.models import Count
import random

class Splitter(AbstractSplitter):
    def __init__(self, namingStub):
        super(Splitter, self).__init__(namingStub)
        self.TEST_PERCENT = 0.30
        self.MAX_PARTITION_SIZE = 35 # Magic # TODO: Please de-magic this.
        self.MIN_TRAIN_SIZE = 10 # WEKA requires at least 10 training points for SVMs.
        self.num_splits = 1

    def split(self, reactions, verbose=False):
        splits = [split(reactions, verbose=verbose) for i in xrange(num_splits)]
        return splits

    def single_split(reactions, verbose=False):
        # Partition the reactions based on what compounds they contain.
        key_counts = self._count_compound_sets(reactions).items()
        random.shuffle(key_counts)

        max_test_size = self.TEST_PERCENT*sum(count for key, count in key_counts)
        test_size = 0
        train_size = 0

        # Determine which partitions should be tested.
        test_keys = []
        for key, count in key_counts:
            if train_size < self.MIN_TRAIN_SIZE:
                train_size += count
            elif test_size < max_test_size and count <= self.MAX_PARTITION_SIZE:
                test_keys.append((key, count))
                test_size += count

        # Get the primary keys of the reactions that should be tested.
        test_reactions = set()
        for key, count in test_keys:
            # Get the reactions that have only the compounds in a given key.
            partition = reactions.annotate(c=Count('compounds')).filter(c=len(key))
            for compound in key:
                partition = partition.filter(compounds=compound)

            test_reactions.update(reaction.pk for reaction in partition)

        # Actually perform the query.
        train = reactions.exclude(pk__in=test_reactions)
        test = reactions.filter(pk__in=test_reactions)

        if verbose:
            print "Split into train ({}), test ({})".format(train.count(), test.count())

        return (self.package(train), self.package(test))

    def _count_compound_sets(self, reactions):
        compound_sets = [frozenset(Compound.objects.filter(reaction=rxn)) for rxn in reactions]
        counts = {b: sum(a==b for a in compound_sets) for b in compound_sets}
        return counts

