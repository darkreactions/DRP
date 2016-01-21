#!/usr/bin/env python
'''A module containing tests for the Balanced Classification Rate SVM'''

import unittest
from DRP.models import PerformedReaction, ModelContainer, Descriptor
from decorators import createsPerformedReactionSet
from DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsPerformedReactionSet
class BCRWekaSVM(DRPTestCase):

  def runTest(self):

    reactions = PerformedReaction.objects.all()

    container = ModelContainer("weka_BCR", "SVM_BCR", splitter="KFoldSplitter",
                               reactions=reactions)
    container.save()

    predictors = Descriptor.objects.filter(heading="testNumber")
    responses = Descriptor.objects.filter(heading="outcome")
    container.build(predictors, responses)

    #TODO: We should test the ModelContainer "predict" method here as well.

    print container.summarize()

suite = unittest.TestSuite([
          loadTests(BCRWekaSVM),
          ])

if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
