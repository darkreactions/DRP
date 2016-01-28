#!/usr/bin/env python
'''A module containing tests for the J48 implementation of the C4.5 classifier'''

import unittest
from DRP.models import PerformedReaction, ModelContainer, Descriptor
from decorators import createsPerformedReactionSetOrd, createsPerformedReactionSetBool
from DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase


class J48(DRPTestCase):

  def runTest(self):

    reactions = PerformedReaction.objects.all()

    container = ModelContainer("weka", "J48", splitter="KFoldSplitter",
                               reactions=reactions)
    container.save()

    predictors = Descriptor.objects.filter(heading="testNumber")
    responses = Descriptor.objects.filter(heading="outcome")
    container.build(predictors, responses)

    #TODO: We should test the ModelContainer "predict" method here as well.

    #print container.summarize()


J48Bool = createsPerformedReactionSetBool(J48)


suite = unittest.TestSuite([
          loadTests(J48Bool)
          ])


if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
