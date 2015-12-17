#!/usr/bin/env python
'''A module containing tests for the KFoldSplitter '''

import unittest
from decorators import createsPerformedReactionSet
from DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase
from DRP.models import Reaction
from DRP.ml_models.splitters.KFoldSplitter import Splitter

@createsPerformedReactionSet
class BasicWekaSVM(DRPTestCase):

  def runTest(self):
    splitterObj = Splitter()

    reactions = Reaction.objects.all()

    splits = splitterObj.split(reactions)

    # Make sure we have K different cross-fold splits.
    assert( splitterObj.k == len(splits) )

    # Make sure each split is not empty.
    for train, test in splits:
      assert( len(train) > 0 )
      assert( len(test) > 0 )


suite = unittest.TestSuite([
          loadTests(BasicWekaSVM),
          ])

if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
