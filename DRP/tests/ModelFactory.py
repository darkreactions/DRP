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
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH"], ["Org"], [0.13],
                          {"outcome":1, "testNumber":5.04})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH"], ["Org"], [0.71],
                          {"outcome":1, "testNumber":5.0})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH"], ["Org"], [0.18],
                          {"outcome":1, "testNumber":5.0})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH"], ["Org"], [0.23],
                          {"outcome":1, "testNumber":6.0})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH"], ["Org"], [0.14],
                          {"outcome":1, "testNumber":5.0})
@createsPerformedReaction("Watchmen", "Rorschach",["dmed"], ["Org"], [0.1],
                          {"outcome":1, "testNumber":5.0})
@createsPerformedReaction("Watchmen", "Rorschach",["dmed"], ["Org"], [0.71],
                          {"outcome":1, "testNumber":5.08})
@createsPerformedReaction("Watchmen", "Rorschach",["dmed"], ["Org"], [0.18],
                          {"outcome":1, "testNumber":5.0})
@createsPerformedReaction("Watchmen", "Rorschach",["dmed"], ["Org"], [0.23],
                          {"outcome":1, "testNumber":6.0})
@createsPerformedReaction("Watchmen", "Rorschach",["dmed"], ["Org"], [0.14],
                          {"outcome":1, "testNumber":5.07})
@createsPerformedReaction("Watchmen", "Rorschach",["dabco"], ["Org"], [0.1],
                          {"outcome":1, "testNumber":5.06})
@createsPerformedReaction("Watchmen", "Rorschach",["dabco"], ["Org"], [0.13],
                          {"outcome":1, "testNumber":5.05})
@createsPerformedReaction("Watchmen", "Rorschach",["dabco"], ["Org"], [0.71],
                          {"outcome":1, "testNumber":5.04})
@createsPerformedReaction("Watchmen", "Rorschach",["dabco"], ["Org"], [0.23],
                          {"outcome":1, "testNumber":6.02})
@createsPerformedReaction("Watchmen", "Rorschach",["dabco"], ["Org"], [0.14],
                          {"outcome":1, "testNumber":5.01})

# And create some more complex ones...
@createsPerformedReaction("Watchmen", "Rorschach",["Water"], ["Water"], [0.53],
                          {"outcome":1, "testNumber":7.0})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH", "dmed", "Water"],
                          ["Org", "Org", "Water"],[0.31, 0.3, 0.5],
                          {"outcome":4, "testNumber":0.1})
@createsPerformedReaction("Watchmen", "Rorschach",["dabco", "Water"],
                          ["Org", "Water"],[0.3, 0.32],
                          {"outcome":3, "testNumber":0.1})
@createsPerformedReaction("Watchmen", "Rorschach",["dmed", "Water"],
                          ["Org", "Water"],[0.3, 0.34],
                          {"outcome":3, "testNumber":0.1})
@createsPerformedReaction("Watchmen", "Rorschach",["EtOH", "Water"],
                          ["Org", "Water"],[0.3, 0.35],
                          {"outcome":3, "testNumber":0.1})
@createsPerformedReaction("Watchmen", "Rorschach",["dmed", "Water"],
                          ["Org", "Water"],[0.4, 0.56],
                          {"outcome":4, "testNumber":0.01})
@createsPerformedReaction("Watchmen", "Rorschach",["Water", "EtOH"],
                          ["Water", "Org"],[0.4, 0.57],
                          {"outcome":4, "testNumber":0.02})
class BasicWekaSVM(DRPTestCase):

  def runTest(self):
    reactions = PerformedReaction.objects.all()

    predictors = ["testNumber"]
    responses = ["outcome"]

    factory = ModelFactory()
    factory.build(reactions, predictors, responses,
                  modelLibrary="weka", modelType="svm",
                  debug=True)


suite = unittest.TestSuite([
          loadTests(BasicWekaSVM),
          ])

if __name__=='__main__':
  runTests(suite)
  #Runs the test- a good way to check that this particular test set works without having to run all the tests.
