"""A splitter to create training and test sets dependent upon apparently disctinct chemistry."""
from .abstractSplitter import AbstractSplitter
from DRP.models.compound import Compound
from DRP.models import CatRxnDescriptor
from django.db.models import Count
import random
import logging
logger = logging.getLogger(__name__)


class Splitter(AbstractSplitter):
    """
    The splitter class.

    This uses the reaction hash descriptor calculated from the reactant names to determine
    whether or not reactions have distinct chemistry.
    """

    def __init__(self, namingStub, num_splits=1, margin_percent=0.01, test_percent=0.33):
        """Specify the number of splits, naming pattern, and ratio of test/traning dataset sizes."""
        super(Splitter, self).__init__(namingStub)
        self.test_percent = test_percent
        self.margin_percent = margin_percent
        self.num_splits = num_splits

    def split(self, reactions, verbose=False):
        """Actually perform the split."""
        super(Splitter, self).split(reactions, verbose=verbose)
        key_counts = self._count_compound_sets(reactions).items()
        splits = [self._single_split(
            reactions, key_counts, verbose=verbose) for i in xrange(self.num_splits)]
        return splits

    def _single_split(self, reactions, key_counts, verbose=False):
        random.shuffle(key_counts)

        total_size = sum(count for key, count in key_counts)
        goal_size = self.test_percent * total_size
        margin = total_size * self.margin_percent
        test_size = 0

        # Determine which partitions should be tested.
        test_keys = []
        for key, count in key_counts:
            if (abs(test_size + count - goal_size) < abs(test_size - goal_size)) and (test_size + count < goal_size + margin):
                # Adding the next set gets us closer to the goal size without
                # going over the error margin. Do it
                test_size += count
                test_keys.append(key)
            elif abs(test_size - goal_size) < margin:
                # We shouldn't add the next set and we're within the error
                # margin. We're done
                break
            # The next set is too big and we're not within the error margin.
            # Try skipping it
        else:
            raise RuntimeError(
                'Failed to make a split under the given parameters.')

        rxnhash_descriptor = CatRxnDescriptor.objects.get(
            heading='rxnSpaceHash1')
        rxn_hashes = reactions.filter(
            catrxndescriptorvalue__descriptor=rxnhash_descriptor)
        test = rxn_hashes.filter(catrxndescriptorvalue__value__in=test_keys)
        train = rxn_hashes.exclude(catrxndescriptorvalue__value__in=test_keys)

        if verbose:
            logger.info("Split into train ({}), test ({})".format(
                train.count(), test.count()))

        return (self.package(train), self.package(test))

    def _count_compound_sets(self, reactions):
        rxnhash_descriptor = CatRxnDescriptor.objects.get(
            heading='rxnSpaceHash1')
        rxn_hashes = reactions.filter(
            catrxndescriptorvalue__descriptor=rxnhash_descriptor)
        counts = {rxnhash_val: rxn_hashes.filter(catrxndescriptorvalue__value=rxnhash_val).count(
        ) for rxnhash_val in rxn_hashes.values_list('catrxndescriptorvalue__value', flat=True).distinct()}
        return counts
