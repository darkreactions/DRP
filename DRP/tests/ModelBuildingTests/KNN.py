#!/usr/bin/env python
'''A module containing tests for the K nearest neighbors classifier'''

import unittest
from DRP.models import PerformedReaction, ModelContainer, Descriptor
from DRP.tests.decorators import createsPerformedReactionSetOrd, createsPerformedReactionSetBool
from DRP.tests.DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase


class KNN(DRPTestCase):

  def runTest(self):

    reactions = PerformedReaction.objects.all()

    container = ModelContainer("weka", "KNN", splitter="KFoldSplitter",
                               reactions=reactions)
    container.save()

    predictors = Descriptor.objects.filter(heading="testNumber")
    responses = Descriptor.objects.filter(heading="outcome")
    container.build(predictors, responses)

    #TODO: We should test the ModelContainer "predict" method here as well.

    #print container.summarize()


KNNBool = createsPerformedReactionSetBool(KNN)


suite = unittest.TestSuite([
          loadTests(KNNBool)
          ])


if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
