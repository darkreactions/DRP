#!/usr/bin/env python
'''A module containing tests for the ModelFactory class'''

import unittest
from DRP.models import PerformedReaction
from DRP.ml_models.ModelFactory import ModelFactory
from decorators import createsPerformedReaction, createsCompound, joinsLabGroup, createsChemicalClass, createsUser, createsCompoundRole, createsRxnDescriptor
from DRPTestCase import DRPTestCase, runTests
loadTests = unittest.TestLoader().loadTestsFromTestCase

@createsUser('Rorschach', 'whatareyouwaitingfor')
@joinsLabGroup('Rorschach', 'Watchmen')
@createsChemicalClass('Org', 'Organic')
@createsChemicalClass('Water', 'Water')
@createsCompound('EtOH', 682, 'Org', 'Watchmen')
@createsCompound('dmed', 67600, 'Org', 'Watchmen')
@createsCompound('dabco', 8882, 'Org', 'Watchmen')
@createsCompound('Water', 937, 'Water', 'Watchmen')
@createsCompoundRole('Org', 'Organic')
@createsCompoundRole('Water', 'Water')
@createsRxnDescriptor("outcome", "OrdRxnDescriptor", options={"maximum":4, "minimum":1})
@createsRxnDescriptor("testNumber", "NumRxnDescriptor")

# Create a bunch of simple sample reactions.
@createsPerformedReaction("r1", "Watchmen", "Rorschach",["EtOH"], ["Org"], [0.13],
                          {"outcome":1, "testNumber":5.04})
@createsPerformedReaction("r2", "Watchmen", "Rorschach",["EtOH"], ["Org"], [0.71],
                          {"outcome":1, "testNumber":5.0})
@createsPerformedReaction("r3", "Watchmen", "Rorschach",["dmed"], ["Org"], [0.1],
                          {"outcome":1, "testNumber":5.0})
@createsPerformedReaction("r4", "Watchmen", "Rorschach",["dmed"], ["Org"], [0.18],
                          {"outcome":1, "testNumber":5.0})
@createsPerformedReaction("r5", "Watchmen", "Rorschach",["dabco"], ["Org"], [0.71],
                          {"outcome":1, "testNumber":5.04})
@createsPerformedReaction("r6", "Watchmen", "Rorschach",["dabco"], ["Org"], [0.14],
                          {"outcome":1, "testNumber":5.01})
@createsPerformedReaction("r7", "Watchmen", "Rorschach",["Water"], ["Water"], [0.14],
                          {"outcome":1, "testNumber":5.01})
@createsPerformedReaction("r8", "Watchmen", "Rorschach",["EtOH", "dmed", "Water"],
                          ["Org", "Org", "Water"],[0.31, 0.3, 0.5],
                          {"outcome":4, "testNumber":0.1})
@createsPerformedReaction("r9", "Watchmen", "Rorschach",["dmed", "EtOH", "Water"],
                          ["Org", "Org", "Water"],[0.2, 0.34, 0.3],
                          {"outcome":4, "testNumber":0.1})
@createsPerformedReaction("r10", "Watchmen", "Rorschach",["dabco", "Water"],
                          ["Org", "Water"],[0.3, 0.32],
                          {"outcome":3, "testNumber":0.1})
@createsPerformedReaction("r11", "Watchmen", "Rorschach",["Water", "dabco"],
                          ["Water", "Org"],[0.3, 0.34],
                          {"outcome":3, "testNumber":0.1})
@createsPerformedReaction("r12", "Watchmen", "Rorschach",["dmed", "Water"],
                          ["Org", "Water"],[0.3, 0.34],
                          {"outcome":3, "testNumber":0.1})
@createsPerformedReaction("r13", "Watchmen", "Rorschach",["EtOH", "Water"],
                          ["Org", "Water"],[0.3, 0.35],
                          {"outcome":3, "testNumber":0.1})
@createsPerformedReaction("r14", "Watchmen", "Rorschach",["dmed", "Water"],
                          ["Org", "Water"],[0.4, 0.56],
                          {"outcome":4, "testNumber":0.01})
@createsPerformedReaction("r15", "Watchmen", "Rorschach",["Water", "EtOH"],
                          ["Water", "Org"],[0.4, 0.57],
                          {"outcome":4, "testNumber":0.02})
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
