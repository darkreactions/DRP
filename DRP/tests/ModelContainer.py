#!/usr/bin/env python
'''A module containing tests for the ModelFactory class'''

import unittest
from DRP.models import PerformedReaction
from decorators import createsPerformedReactionSet
from DRP.models import PerformedReaction, ModelContainer, Descriptor
from decorators import createsPerformedReaction, createsCompound, joinsLabGroup, createsChemicalClass, createsUser, createsCompoundRole, createsRxnDescriptor
from DRPTestCase import DRPTestCase, runTests
from DRP.ml_models.splitters.KFoldSplitter import Splitter
loadTests = unittest.TestLoader().loadTestsFromTestCase
 
@createsUser('Rorschach', 'whatareyouwaitingfor')
@joinsLabGroup('Rorschach', 'Watchmen')
@signsExampleLicense("Rorschach")
@createsPerformedReactionSet
class BasicWekaSVM(DRPTestCase):

  def runTest(self):
    predictors = Descriptor.objects.filter(heading="testNumber")
    responses = Descriptor.objects.filter(heading="outcome")

    container = ModelContainer(library="weka", tool="svm",
                               splitter="KFoldSplitter")

    reactions = PerformedReaction.objects.all()

    splitter = Splitter()
    for training, test in splitter.split(reactions):
      container.build(training, test, predictors, responses)

    container.save()
    print container.summarize()


suite = unittest.TestSuite([
          loadTests(BasicWekaSVM),
          ])

if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
