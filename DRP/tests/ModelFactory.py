#!/usr/bin/env python
'''A module containing tests for the ModelFactory class'''

import unittest
from DRP.models import PerformedReaction
from DRP.ml_models.ModelFactory import ModelFactory
from decorators import createsPerformedReactionSet
from DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase
 
@createsUser('Rorschach', 'whatareyouwaitingfor')
@joinsLabGroup('Rorschach', 'Watchmen')
@signsExampleLicense("Rorschach")
@createsPerformedReactionSet
class BasicWekaSVM(DRPTestCase):

  def runTest(self):
    reactions = PerformedReaction.objects.all()

    predictors = ["testNumber"]
    responses = ["outcome"]

    factory = ModelFactory()
    model = factory.build(reactions, predictors, responses,
                          modelLibrary="weka", modelType="svm",
                          debug=True)
    model.summarize()


suite = unittest.TestSuite([
          loadTests(BasicWekaSVM),
          ])

if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
