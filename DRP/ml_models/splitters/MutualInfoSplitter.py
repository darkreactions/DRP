from AbstractSplitter import AbstractSplitter
from DRP.models.Compound import Compound
from django.db.models import Count
import random
import warnings

class Splitter(AbstractSplitter):
    def __init__(self, namingStub, num_splits=1, max_partition_size=35, test_percent=0.3, min_train_size=10):
        super(Splitter, self).__init__(namingStub)
        self.test_percent = test_percent
        self.max_partition_size = max_partition_size
        self.min_train_size = min_train_size
        self.num_splits = num_splits

        if min_train_size < 10:
            warnings.warn('Min train size is only {}. Weka requires at least 10 training points for SVMs')

    def split(self, reactions, verbose=False):
        key_counts = self._count_compound_sets(reactions).items()
        splits = [self._single_split(reactions, key_counts, verbose=verbose) for i in xrange(self.num_splits)]
        return splits

    def _single_split(self, reactions, key_counts, verbose=False):
        # Partition the reactions based on what compounds they contain.
        random.shuffle(key_counts)

        max_test_size = self.test_percent*sum(count for key, count in key_counts)
        test_size = 0
        train_size = 0

        # Determine which partitions should be tested.
        test_keys = []
        for key, count in key_counts:
            if train_size < self.min_train_size:
                train_size += count
            elif test_size < max_test_size and count <= self.max_partition_size:
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

